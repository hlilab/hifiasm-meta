#include <stdio.h>
#include <stdlib.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Overlaps.h"
#include "ksort.h"
#include "Process_Read.h"
#include "CommandLines.h"
#include "Hash_Table.h"
#include "Correct.h"
#include "Purge_Dups.h"
#include <assert.h>
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>  // debug, for printing uint64
#include "Overlaps_hamt.h"
#include "kthread.h"

KRADIX_SORT_INIT(hamtu64, uint64_t, uint64_t, 8)  // hamt

// (hifiasm-meta: when to sort and use binary search over linear search)
// note: sorted/not sorted is indicated by global variables under asm_opt.
//       Must use `add_ma_hit_t_alloc` to add entries to buffer.
#define LINEAR_BF_BREAKEVEN_POINT 64
void sort_paf_buffers_by_targetID(ma_hit_t_alloc *source, long long n_read);
void sort_one_paf_buffer_by_targetID(ma_hit_t_alloc *x);


// (hifiasm: filter short potential unitigs)
#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3

// (hifiasm)
#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

uint32_t debug_purge_dup = 0;

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)

#define ma_hit_key_tn(a) ((a).tn)
KRADIX_SORT_INIT(hit_tn, ma_hit_t, ma_hit_key_tn, member_size(ma_hit_t, tn))

#define ma_hit_key_qns(a) ((a).qns)
KRADIX_SORT_INIT(hit_qns, ma_hit_t, ma_hit_key_qns, member_size(ma_hit_t, qns))


#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(arch64, uint64_t, generic_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(arch32, uint32_t, generic_key, 4)

///#define Hap_Align_key(a) ((((uint64_t)((a).is_color))<<63)|(((uint64_t)((a).t_id))<<32)|((uint64_t)((a).q_pos)))
///#define Hap_Align_key(a) ((((uint64_t)((a).is_color))<<32)|(((uint64_t)((a).t_id))<<33)|((uint64_t)((a).q_pos)))
#define Hap_Align_key(a) ((((uint64_t)((a).t_id))<<33)|((uint64_t)((a).q_pos)))
KRADIX_SORT_INIT(Hap_Align_sort, Hap_Align, Hap_Align_key, 8)


KSORT_INIT_GENERIC(uint32_t)

///this value has been updated at the first line of build_string_graph_without_clean
long long min_thres;

void ma_hit_sort_tn(ma_hit_t *a, long long n)
{
	radix_sort_hit_tn(a, a + n);
}

void ma_hit_sort_qns(ma_hit_t *a, long long n)
{
	radix_sort_hit_qns(a, a + n);
}

void sort_kvec_t_u64_warp(kvec_t_u64_warp* u_vecs, uint32_t is_descend)
{
    radix_sort_arch64(u_vecs->a.a, u_vecs->a.a + u_vecs->a.n);
    if(is_descend)
    {
        uint64_t i, uInfor;
        for (i = 0; i < (u_vecs->a.n>>1); ++i) 
        {
            uInfor = u_vecs->a.a[i];
            u_vecs->a.a[i] = u_vecs->a.a[u_vecs->a.n - i - 1];
            u_vecs->a.a[u_vecs->a.n - i - 1] = uInfor;
        }
    }
} 



///////////////////////////////////////////////
//  assembly graph / arc / overlap interface //
///////////////////////////////////////////////

asg_t *asg_init(void)
{
	return (asg_t*)calloc(1, sizeof(asg_t));
}

void asg_destroy(asg_t *g)
{
	if (g == 0) return;
	free(g->seq); free(g->idx); free(g->arc); free(g->seq_vis);
     
    if(g->n_F_seq > 0 && g->F_seq)
    {
        uint32_t i = 0;
        for (i = 0; i < g->n_F_seq; i++)
        {
            if(g->F_seq[i].a) free(g->F_seq[i].a);
            if(g->F_seq[i].s) free(g->F_seq[i].s);
        }

        free(g->F_seq);
    }
    

    free(g);
}

void asg_arc_sort(asg_t *g)
{
	radix_sort_asg(g->arc, g->arc + g->n_arc);
}


void add_overlaps(ma_hit_t_alloc* source_paf, ma_hit_t_alloc* dest_paf, uint64_t* source_index, long long listLen)
{
    long long i;
    ma_hit_t* tmp;
    for (i = 0; i < listLen; i++)
    {
        tmp = &(source_paf->buffer[(uint32_t)(source_index[i])]);
        add_ma_hit_t_alloc(dest_paf, tmp);
    }
}


void remove_overlaps(ma_hit_t_alloc* source_paf, uint64_t* source_index, long long listLen)
{
    long long i, m;
    for (i = 0; i < listLen; i++)
    {
        source_paf->buffer[(uint32_t)(source_index[i])].qns = (uint64_t)(-1);
    }

    m = 0;
    for (i = 0; i < source_paf->length; i++)
    {
        if(source_paf->buffer[i].qns != (uint64_t)(-1))
        {
            source_paf->buffer[m] = source_paf->buffer[i];
            m++;
        }
    }
    source_paf->length = m;
}

void add_overlaps_from_different_sources(ma_hit_t_alloc* source_paf_list, ma_hit_t_alloc* dest_paf, 
uint64_t* source_index, long long listLen)
{
    long long i;
    ma_hit_t ele;
    ma_hit_t* tmp;
    uint32_t source_n, source_i;
    for (i = 0; i < listLen; i++)
    {
        source_n = source_index[i] >> 32;
        source_i = (uint32_t)(source_index[i]);
        tmp = &(source_paf_list[source_n].buffer[source_i]);

        ele.del = 0;
        ele.rev = tmp->rev;
        ele.qns = Get_tn((*tmp));
        ele.qns = ele.qns << 32;
        ele.qns = ele.qns | (uint64_t)(Get_ts((*tmp)));
        ele.qe = Get_te((*tmp));

        ele.tn = Get_qn((*tmp));
        ele.ts = Get_qs((*tmp));
        ele.te = Get_qe((*tmp));

        ele.bl = R_INF.read_length[ele.tn];
        ele.ml = tmp->ml;
        ele.el = tmp->el;
        ele.no_l_indel = tmp->no_l_indel;
        
        add_ma_hit_t_alloc(dest_paf, &ele);
    }
}

void ma_ug_destroy(ma_ug_t *ug)
{
	uint32_t i;
	if (ug == 0) return;
	for (i = 0; i < ug->u.n; ++i) {
		free(ug->u.a[i].a);
		free(ug->u.a[i].s);
	}
	free(ug->u.a);
	asg_destroy(ug->g);
	free(ug);
}

uint64_t *asg_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);


    /**
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
	**/
    ///so if we use high 32-bit, we store the index of each qn with two direction
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || a[i-1].ul>>32 != a[i].ul>>32)
			idx[a[i-1].ul>>32] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void asg_arc_index(asg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = asg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

void asg_seq_set(asg_t *g, int sid, int len, int del)
{
	///just malloc size
	if (sid >= (int)g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
	}

	if (sid >= g->n_seq) g->n_seq = sid + 1;
	
    g->seq[sid].del = !!del;
    g->seq[sid].len = len;
}

ma_utg_t* asg_F_seq_set(asg_t *g, int iid)
{
    if(iid < (int)g->r_seq) return NULL;
    uint32_t index = iid - g->r_seq, pre_n_F_seq = g->n_F_seq;

	///just malloc size
	if (index >= g->n_F_seq) {
		g->n_F_seq = index + 1;
		kv_roundup32(g->n_F_seq);
        g->F_seq = (ma_utg_t*)realloc(g->F_seq, g->n_F_seq*sizeof(ma_utg_t));
        memset(g->F_seq + pre_n_F_seq, 0, sizeof(ma_utg_t)*(g->n_F_seq-pre_n_F_seq));
	}

    return &(g->F_seq[index]);
}

// hard remove arcs marked as "del"
void asg_arc_rm(asg_t *g)
{
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
	p->v : |___________31___________|__________1___________|
				        tns              relative strand between query and target
	p->ol: overlap length
	**/
	uint32_t e, n;
	///just clean arc requiring: 1. arc it self must be available 2. both the query and target are available
	for (e = n = 0; e < g->n_arc; ++e) {
		//u and v is the read id
		uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

void asg_cleanup(asg_t *g)
{
    ///remove overlaps, instead of reads
    ///remove edges with del, and free idx
	asg_arc_rm(g);
	if (!g->is_srt) {
		/**
		 * sort by ul, that is, sort by qns + direction
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
		**/
		asg_arc_sort(g);
		g->is_srt = 1;
	}
	///index the overlaps in graph with query id
	if (g->idx == 0) asg_arc_index(g);
}

// delete multi-arcs
/**
 * remove edges like:   v has two out-edges to w
**/
int asg_arc_del_multi(asg_t *g)
{
	//the number of nodes are number of read times 2
	uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
	for (v = 0; v < n_vtx; ++v) {
		///out-nodes of v
		asg_arc_t *av = asg_arc_a(g, v);
		int32_t i, nv = asg_arc_n(g, v);
		///if v just have one out-node, there is no muti-edge
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].v] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (n_multi) asg_cleanup(g);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
    }
	
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int asg_arc_del_asymm(asg_t *g)
{
	uint32_t e, n_asymm = 0;
	///g->n_arc is the number of overlaps 
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
		uint32_t i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].v == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (n_asymm) asg_cleanup(g);
    if(VERBOSE >= 1)
    {
	    fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
    }
	return n_asymm;
}


void asg_symm(asg_t *g)
{
	asg_arc_del_multi(g);
	asg_arc_del_asymm(g);
	g->is_symm = 1;
}

void init_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->size = 0;
    x->buffer = NULL;
    x->length = 0;
}

void clear_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->length = 0;
}

void resize_ma_hit_t_alloc(ma_hit_t_alloc* x, uint32_t size)
{
	if (size > x->size) {
		x->size = size;
		kroundup32(x->size);
		REALLOC(x->buffer, x->size);
	}
}

void destory_ma_hit_t_alloc(ma_hit_t_alloc* x, long long total_reads)
{
	// free(x->buffer);
    for (long long i=0; i<total_reads; i++){
        free(x[i].buffer);
    }
    free(x);
}

void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element)
{
	if (x->length + 1 > x->size) {
		x->size = x->length + 1;
		kroundup32(x->size);
		REALLOC(x->buffer, x->size);
	}
    int qn = (int)(element->qns>>32);
    int tn = (int)(element->tn);
	x->buffer[x->length++] = *element;
}

long long get_specific_overlap_binarysearch(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn){
	// NOTE
    //     Internal of get_specific_overlap, do not use directly.
    
    if (x->length==0) return -1;
    
    // the regulary binary search
	int low=0, high=x->buffer_sorted_by_tn, mid;
    int dbg_cnt = 0;
    while (low<high){
        mid = (high+low)/2;
        assert ((uint32_t)(x->buffer[mid].qns>>32) == qn);  // debug - should be guaranteed but why does ha to check this?
        if (tn<=x->buffer[mid].tn){
            high = mid;
        }else{
            low = mid+1;
        }
        dbg_cnt++;
        if (dbg_cnt>10000){
            fprintf(stderr, "[E::%s] probably infinite loop; array traceback:\n", __func__);
            for (int i=0; i<x->buffer_sorted_by_tn; i++){
                fprintf(stderr, "%d ", (int)x->buffer[i].tn);
            }
            fprintf(stderr, "\n");
            fflush(stderr);
            exit(1);
        }
    }
    if (x->buffer[high].tn!=tn) {
        if (x->length-1>x->buffer_sorted_by_tn){  // linear search the unsorted portion, if any
            for (int i=x->buffer_sorted_by_tn+1; i<x->length; i++){
                if (x->buffer[i].tn==tn)  return i;
            }
        }
        return -1;
    }else{
        return high;
    }
}
long long get_specific_overlap_linearsearch(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn)
{
    uint32_t i;
    for (i = 0; i < x->length; i++){
        if(x->buffer[i].tn == tn && ((uint32_t)(x->buffer[i].qns>>32)) == qn){
            return i; 
        }
    }
    return -1;  // not found
}
long long get_specific_overlap(ma_hit_t_alloc* x, uint32_t qn, uint32_t tn)
{
	if ( (!asm_opt.get_specific_overlap_is_use_bf) || x->length<LINEAR_BF_BREAKEVEN_POINT){
		return get_specific_overlap_linearsearch(x, qn, tn);
	}else{
		return get_specific_overlap_binarysearch(x, qn, tn);
	}
}

ma_hit_t *get_specific_overlap_handle(ma_hit_t_alloc *sources, uint32_t qn, uint32_t tn){  // note: no direction bit
    for (int i=0; i<sources[qn].length; i++){
        // if (sources[qn].buffer[i].del){continue;}
        if (sources[qn].buffer[i].tn == tn){
            return &sources[qn].buffer[i];
        }
    }
    return NULL;
}


int hamt_count_mismatches_in_an_exact_overlap(ma_hit_t *h, UC_Read *read1, UC_Read *read2){
    // FUNC
    //     Given a handle, get the corrected bases and compare.
    // NOTE
    //     Since we can't know indels, this function aborts if there's consecutive mismatches.
    // RETURN
    //     -1 if overlap is not exact, or has long indel, or might have any indel (a rough check)
    //     INT>=0 number of mismatches otherwise
    if (!h->el) {return -1;}
    if (!h->no_l_indel) {return -1;}

    int nm = 0;
    uint32_t v, w, qs;
    v = (uint32_t) (h->qns>>32);
    w = h->tn;
    qs = (uint32_t)h->qns;
    recover_UC_Read(read1, &R_INF, v);
    recover_UC_Read(read2, &R_INF, w);
    uint32_t i1, i2;
    char a, b;
    int prv_is_nm = 0;

    for (uint32_t i=0; i<h->qe-qs; i++){
        if (h->rev){
            i1 = qs+i;
            i2 = h->te-i;
            a = seq_nt4_table[(uint8_t)read1->seq[i1]];
            b = seq_nt4_table[(uint8_t)read1->seq[i2]];
            if (a!=b && a+b==4){  // rev comp
                prv_is_nm = 0;
                continue;
            }else{
                if (prv_is_nm) {return -2;}  // consecutive mismatch
                nm++;
                prv_is_nm = 1;
            }
        }else{
            i1 = qs+i;
            i2 = h->ts+i;
            if (read1->seq[i1]!=read2->seq[i2]){
                if (prv_is_nm) {return -2;}  // consecutive mismatch
                nm++;
                prv_is_nm = 1;
            }else{
                prv_is_nm = 0;
            }
        }
    }

    return nm;
}
void hamt_count_mismatches_in_an_exact_overlap_debug(int ID1, int ID2){
    UC_Read r1, r2;
    ma_hit_t *h;
    uint32_t qs;
    init_UC_Read(&r1);
    init_UC_Read(&r2);
    for (int i=0; i<R_INF.paf[ID1].length; i++){
        h = &R_INF.paf[ID1].buffer[i];
        if (h->tn!=ID2) continue;
        qs = (uint32_t) (h->qns>>32);
        fprintf(stderr, "is rev: %d\n", (int)h->rev);
        fprintf(stderr, "no long indel: %d\n", (int)h->no_l_indel);

        recover_UC_Read(&r1, &R_INF, ID1);
        recover_UC_Read(&r2, &R_INF, ID2);
        for (int idx=0; idx<h->qe-qs; idx++){
            fprintf(stderr, "  %d\t%c\t%c\n", idx, r1.seq[qs+idx], r2.seq[h->te-idx]);
        }

    }
    destory_UC_Read(&r1);
    destory_UC_Read(&r2);

}


int hit_is_perfect_containment(ma_hit_t *hit){
    uint32_t qs = (uint32_t)hit->qns;
    uint32_t qn = (uint32_t)(hit->qns>>32);
    if ( (qs==0 && hit->qe==R_INF.read_length[qn]) ||
         (hit->ts==0 && hit->te==R_INF.read_length[hit->tn]) ){
        return 1;
    }
    return 0;
}
void hamt_try_rescue_containment(ma_hit_t_alloc* sources, long long nb_reads){
    // FUNC
    //     Set reverse overlap for one-direction overlaps that happen to be perfect containments.
    //     The added overlap will be in the same collection as the query.
    // NOTE
    //     Introduced to resolve a problem happens in very high coverage occasions (roughly >200x).
    //     Doesn't align the other direction for double check.
    //     This was done inside `normalize_ma_hit_t_single_side_advance`. Seperating out
    //       for clarity.
    //     The better way is to either a) verify this rescued overlap, or b) probe from the 
    //       final ov. b) was clumsy. Do we need a)? (r27)
    uint32_t qID, tID;
    int rev_index;
    ma_hit_t new_hit;
    for (int i=0; i<nb_reads; i++){
        for (int j=0; j<sources[i].length; j++){
            qID = Get_qn(sources[i].buffer[j]);
            tID = Get_tn(sources[i].buffer[j]);
            rev_index = get_specific_overlap(&(sources[tID]), tID, qID);  // find ovlp of t->q
            if ( (rev_index==-1) && (hit_is_perfect_containment(&sources[i].buffer[j])) ){
                set_reverse_overlap(&new_hit, &(sources[i].buffer[j]));  // add the reversed overlap
                add_ma_hit_t_alloc(&(sources[tID]), &new_hit);
                if (asm_opt.get_specific_overlap_is_use_bf){
                    sort_one_paf_buffer_by_targetID(&sources[tID]);
                }
            }
        }
    }
}

void hamt_smash_haplotype(ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, long long nb_reads){
    // FUNC
    //    Copy het overlaps into the hom buffer. Will not diallocate reverse_sources.
    // NOTE
    //    Introduced after tests in the marine sample. If the input has very low coverage (<10x everywhere),
    //       and also contains many haplotypes, maybe ignore phasing for most of the time
    //       can help to build a better graph.
    uint32_t qID, tID;
    int index;
    ma_hit_t *hit;
    for (int i=0; i<nb_reads; i++){
        for (int j=0; j<reverse_sources[i].length; j++){
            hit = &reverse_sources[i].buffer[j];
            qID = Get_qn(*hit);
            tID = Get_tn(*hit);

            // check if the overlap is already present in the `sources` (e.g. added by a resuce)
            index = get_specific_overlap(&(sources[qID]), qID, tID);
            if (index!=-1) continue;

            // copy the overlap
            add_ma_hit_t_alloc(&(sources[qID]), hit);
            if (asm_opt.get_specific_overlap_is_use_bf){
                sort_one_paf_buffer_by_targetID(&sources[qID]);
            }
        }
    }
}

uint8_t hamt_check_if_reads_overlap(uint32_t v, uint32_t w, ma_hit_t_alloc *paf,
                                    ma_sub_t *coverage_cut, int allow_del){
    // FUNC
    //    Check if v and w overlap with each other in the given overlap collection (paf).
    // INPUT
    //    coverage_cut could be NULL; if so, won't check the corresponding .del status.
    // TODO
    //    Allow excluding some cases, e.g. containment (based on hit2arc)
    // RETURN
    //    0 N0
    //    1 asymmetric, only v->w (v be the query and w be the target)
    //    2 asymmetric, only w->v 
    //    3 symmetric, both v->w and w->v are found
    uint8_t ret = 0;
    ma_hit_t_alloc *h = &paf[v];
    ma_hit_t *hit = NULL;
    // direction 1
    for (int i=0; i<h->length; i++){
        hit = &(h->buffer[i]);
        if (!allow_del && hit->del) continue;
        if (coverage_cut){
            if (coverage_cut[Get_qn(*hit)].del || coverage_cut[Get_tn(*hit)].del) {
                continue;
            }
        }
        if (hit->tn==w){
            ret |= 1;
        }
    }
    // direction 2
    h = &paf[w];
    for (int i=0; i<h->length; i++){
        hit = &(h->buffer[i]);
        if (!allow_del && hit->del) continue;
        if (coverage_cut){
            if (coverage_cut[Get_qn(*hit)].del || coverage_cut[Get_tn(*hit)].del) {
                continue;
            }
        }
        if (hit->tn==v){
            ret |= 2;
        }
    }
    return ret;
}

// NOT TESTED
int hamt_count_haplotypes_given_reads(ma_hit_t_alloc *reverse_sources,
                                      uint32_t *buf, int buf_len, int allow_del){
    // FUNC
    //     Given an array of read (IDs, no direction bit),
    //       count how many haplotypes are there based on phasing hints of sources/reverse_sources.
    //     Note that this assumes reads in buf are largely overlapping with 
    //       each other - we just peek if there's arcs in reverse paf for positive haplotype hint.
    // RETURN
    //     Number of haplotypes. 0 if buf is empty.
    if (buf_len<=1) {return buf_len;}

    uint32_t *buf_hap = (uint32_t*)malloc(sizeof(uint32_t) * buf_len);
    buf_hap[0] = buf[0];
    int buf_hap_len = 1;
    uint32_t v, w;

    for (int i=1; i<buf_len; i++){
        v = buf[i];  // the incoming candidate
        for (int j=0; j<buf_hap_len; j++){
            w = buf_hap[j];  // from the established collection
            if (hamt_check_if_reads_overlap(v, w, reverse_sources, NULL, allow_del)>0){
                buf_hap[buf_hap_len] = v;
                buf_hap_len++;
                break;
            }
        }
    }

    free(buf_hap);
    return buf_hap_len;
}


int hamt_read_het_to_all_given_reads(ma_hit_t_alloc *sources, uint32_t r,
                                      uint32_t *buf, int buf_len, int allow_del){
    // NOTE
    //     Won't care if r isn't overlapped with reads in the buffer,
    //      in this case the return is positive.
    // RETURN
    //     -1 unknown because empty buffer
    //     0 no
    //     1 yes
    if (buf_len==0) {return -1;}


    for (int i=0; i<buf_len; i++){
        if (hamt_check_if_reads_overlap(buf[i], r, sources, NULL, allow_del)>0){
            return 0;
        }
    }
    return 1;

}

int hamt_read_is_contained_by_others(ma_hit_t_alloc* paf, uint32_t v, uint64_t v_len, 
                                     uint32_t **buf0, int *buf_size0, int *buf_len0){
    // FUNC
    //     Check if v is contained by any other reads in `paf`.
    //     If `buf` is not NULL, store parents in it. 
    // NOTE
    //     For simplicity we only check in the directionv->target here.
    // PAR
    //     v_len: read v's  length
    // RETURN
    //     count of parent reads (0 means not contained)
    int ret = 0;
    ma_hit_t *h;
    uint32_t *buf;
    int buf_size, buf_len;
    if (buf0) {
        buf = *buf0; 
        buf_size = *buf_size0;
        buf_len = *buf_len0;
        assert(buf_size>0);
        assert(buf_len<=buf_size);
    }else{
        buf = 0;
    }
    


    for (int i=0; i<paf[v].length; i++){
        h = &paf[v].buffer[i];
        if (h->del) continue;
        if ((uint32_t)h->qns==0 && h->qe==v_len){  // query is contained by h->tn
            if (buf){  // push the parentID
                if (buf_len==buf_size){
                    buf_size = buf_size + ((buf_size)>>1);
                    buf = (uint32_t*)realloc(buf, buf_size * sizeof(uint32_t) );
                    assert(buf);
                }
                // fprintf(stderr, "len %d, size %d\n", buf_len, buf_size);
                buf[buf_len] = h->tn;
                buf_len = buf_len + 1;
            }
            ret+=1;
        }
    }

    if (buf){
        *buf_size0 = buf_size;
        *buf_len0 = buf_len;
        *buf0 = buf;
    }

    return ret;
}


void normalize_ma_hit_t_single_side_advance(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();
    int verbose = 0;
    int sancheck_type1 = 0;
    int sancheck_type2 = 0;
    int sancheck_type3 = 0;

    long long i, j, index;
    uint32_t qn, tn, is_del = 0;
    long long qLen_0, qLen_1;
    ma_hit_t ele;
    for (i = 0; i < num_sources; i++)  // iterate over vertices
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);

            ///if(sources[i].buffer[j].del) continue;

            index = get_specific_overlap(&(sources[tn]), tn, qn);  // find ovlp of t->q


            ///if(index != -1 && sources[tn].buffer[index].del == 0)
            if(index != -1)  // both t->q and q->t are present
            {
                is_del = 0;
                if(sources[i].buffer[j].del || sources[tn].buffer[index].del)
                {
                    is_del = 1;
                }

                qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);

                if(qLen_0 == qLen_1)
                {
                    ///qn must be not equal to tn
                    ///make sources[qn] = sources[tn] if qn > tn
                    if(qn < tn)  // overwrite the entry ordered higher (later)
                    {
                        set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                        sancheck_type1++;
                    }
                }
                else if(qLen_0 > qLen_1)  // overwrite the reverse entry
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    sancheck_type2++;
                }
                
                sources[i].buffer[j].del = is_del;
                sources[tn].buffer[index].del = is_del;
            }
            else ///means this edge just occurs in one direction
            {  
                set_reverse_overlap(&ele, &(sources[i].buffer[j]));  // add the reversed overlap
                sources[i].buffer[j].del = ele.del = 1;  // but also delete it; note: assuming containment problem in high coverage occasions are already handled
                add_ma_hit_t_alloc(&(sources[tn]), &ele);
                sancheck_type3++;
                if (asm_opt.get_specific_overlap_is_use_bf){
                    sort_one_paf_buffer_by_targetID(&sources[tn]);
                }
            }
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs; 3 types: %d, %d, %d\n\n", __func__, Get_T()-startTime,
                    sancheck_type1, sancheck_type2, sancheck_type3);
    }
}

typedef struct {
    // worker will collect query-target pairs that shall go into `set_reverse_overlap`
    ma_hit_t ***buf1;  // 1st in the pair, array of pointers each points to an array of pointers
    ma_hit_t ***buf2;  // 2nd
    uint32_t **is_hard_add;  // this slot needs a placeholder 1st, it's adding an entry rather than modify the 1st given the 2nd
    uint64_t **is_hard_add_buf;
    // ^ see block with ha comment "///means this edge just occurs in one direction"
    uint64_t *n;  // buffer's current lengths, per thread
    uint64_t *m;  // buffer's current sizes, per thread
    ma_hit_t_alloc* sources;
    long long i_start;
    long long n_read;
}normhit_aux_t;

static void hamt_normalize_ma_hit_t_single_side_advance_worker(void *data, long i_r_, int tid){  // callback for kt_for()
    normhit_aux_t *d = (normhit_aux_t*)data;
    ma_hit_t_alloc *sources = d->sources;
    long i_r = i_r_ + d->i_start;  // batched kt_for
    if (i_r>=d->n_read) return;

    ma_hit_t **buf1 = d->buf1[tid];
    ma_hit_t **buf2 = d->buf2[tid];

    uint32_t *is_hard_add = d->is_hard_add[tid];
    uint64_t *is_hard_add_buf = d->is_hard_add_buf[tid];
    uint64_t *buf_n = &d->n[tid];
    uint64_t *buf_m = &d->m[tid];

    long long i = (long long) i_r; 
    long long j, index;
    uint32_t qn, tn, is_del = 0;
    long long qLen_0, qLen_1;
    ma_hit_t ele;
    for (j = 0; j < sources[i].length; j++)  // iterate of overlaps of this vertex. (.length is Cigar_record's .length)
    {
        if (*buf_n==*buf_m){
            // fprintf(stderr, "expand, tid %d, n %d m%d\n", tid, (int)(*buf_n), (int)(*buf_m));fflush(stderr);
            *buf_m = (*buf_m) + ((*buf_m)>>1);
            d->buf1[tid] = (ma_hit_t**)realloc(d->buf1[tid], (*buf_m)*sizeof(ma_hit_t*));
            d->buf2[tid] = (ma_hit_t**)realloc(d->buf2[tid], (*buf_m)*sizeof(ma_hit_t*));
            d->is_hard_add[tid] = (uint32_t*)realloc(d->is_hard_add[tid], (*buf_m)*sizeof(uint32_t));
            d->is_hard_add_buf[tid] = (uint64_t*)realloc(d->is_hard_add_buf[tid], (*buf_m)*sizeof(uint64_t));
            buf1 = d->buf1[tid];
            buf2 = d->buf2[tid];
            is_hard_add = d->is_hard_add[tid];
            is_hard_add_buf = d->is_hard_add_buf[tid];
        }

        qn = Get_qn(sources[i].buffer[j]);
        tn = Get_tn(sources[i].buffer[j]);

        sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
        index = get_specific_overlap(&(sources[tn]), tn, qn);  // find ovlp of t->q

        ///if(index != -1 && sources[tn].buffer[index].del == 0)
        if(index != -1)  // both t->q and q->t are present
        {
            is_del = 0;
            if(sources[i].buffer[j].del || sources[tn].buffer[index].del)
            {
                is_del = 1;
            }

            qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
            qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);

            if(qLen_0 == qLen_1)
            {
                ///qn must be not equal to tn
                ///make sources[qn] = sources[tn] if qn > tn
                if(qn < tn)  // overwrite the entry ordered higher (later)
                {
                    // set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    buf1[*buf_n] = &(sources[tn].buffer[index]);
                    buf2[*buf_n] = &(sources[i].buffer[j]);
                    is_hard_add[*buf_n] = 0;
                    *buf_n = (*buf_n) + 1;
                }
            }
            else if(qLen_0 > qLen_1)  // overwrite the reverse entry
            {
                // set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                buf1[*buf_n] = &(sources[tn].buffer[index]);
                buf2[*buf_n] = &(sources[i].buffer[j]);
                is_hard_add[*buf_n] = 0;
                *buf_n = (*buf_n) + 1;
            }
            
            sources[i].buffer[j].del = is_del;
            sources[tn].buffer[index].del = is_del;
        }
        else ///means this edge just occurs in one direction
        {  
            // set_reverse_overlap(&ele, &(sources[i].buffer[j]));  // add the reversed overlap
            sources[i].buffer[j].del = 1;  // but also delete it
            // add_ma_hit_t_alloc(&(sources[tn]), &ele);

            buf1[*buf_n] = &(sources[tn].buffer[index]);  // is a placeholder
            buf2[*buf_n] = &(sources[i].buffer[j]);
            is_hard_add[*buf_n] = (tn<<1) | 1; 
            is_hard_add_buf[*buf_n] = (((uint64_t)i)<<32) | ((uint64_t)j);
            *buf_n = (*buf_n) + 1;
        }
    }


}
void hamt_normalize_ma_hit_t_single_side_advance(ma_hit_t_alloc* sources, long long n_read){
    // FUNC
    //     Threaded normalize_ma_hit_t_single_side_advance.
    double startTime = Get_T();
    int verbose = 0;
    int n_cpu = asm_opt.thread_num;
    
    int batch_size = n_cpu*64;
    long long nb_batch =  n_read/((long long)batch_size)+1;
    uint64_t sancheckA = 0, sancheckB = 0;
    

    // init data
    int buf_size = 1024;
    normhit_aux_t aux;
    aux.sources = sources;
    aux.i_start = 0;
    aux.n_read = n_read;
    aux.n = (uint64_t*)calloc(n_cpu, sizeof(uint64_t));
    aux.m = (uint64_t*)calloc(n_cpu, sizeof(uint64_t));
    aux.is_hard_add = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));
    aux.is_hard_add_buf = (uint64_t**)calloc(n_cpu, sizeof(uint64_t*));
    aux.buf1 = (ma_hit_t***)calloc(n_cpu, sizeof(ma_hit_t**));
    aux.buf2 = (ma_hit_t***)calloc(n_cpu, sizeof(ma_hit_t**));
    for (int i=0; i<n_cpu; i++){
        aux.m[i] = buf_size;
        aux.is_hard_add[i] = (uint32_t*)calloc(buf_size, sizeof(uint32_t));
        aux.is_hard_add_buf[i] = (uint64_t*)calloc(buf_size, sizeof(uint64_t));
        aux.buf1[i] = (ma_hit_t**)calloc(buf_size, sizeof(ma_hit_t*));
        aux.buf2[i] = (ma_hit_t**)calloc(buf_size, sizeof(ma_hit_t*));
    }

    fprintf(stderr, "[debug::%s] nb_batch: %d\n", __func__, (int)nb_batch);
    for (long long i_batch=0; i_batch<nb_batch; i_batch++){
        for (int i=0; i<n_cpu; i++){aux.n[i] = 0;}  // reset all buffers
        aux.i_start = i_batch * ((long long)batch_size);

        // collect
        kt_for(n_cpu, hamt_normalize_ma_hit_t_single_side_advance_worker, &aux, (long)batch_size);

        // treat
        ma_hit_t ele;
        for (int tid=0; tid<n_cpu; tid++){
            for (uint64_t n=0; n<aux.n[tid]; n++){
                if (!aux.is_hard_add[tid][n]){
                    set_reverse_overlap(aux.buf1[tid][n], aux.buf2[tid][n]);
                    sancheckA++;
                }
            }
        }
        uint32_t i_tmp, j_tmp;
        for (int tid=0; tid<n_cpu; tid++){
            for (uint64_t n=0; n<aux.n[tid]; n++){
                if (aux.is_hard_add[tid][n]){
                    i_tmp = aux.is_hard_add_buf[tid][n]>>32;
                    j_tmp = (uint32_t)aux.is_hard_add_buf[tid][n];
                    set_reverse_overlap(&ele, &sources[i_tmp].buffer[j_tmp]);
                    ele.del = 1;
                    add_ma_hit_t_alloc(&(sources[aux.is_hard_add[tid][n]>>1]), &ele);
                    sancheckB++;
                }
            }
        }

        // maintain sort order
        if (asm_opt.get_specific_overlap_is_use_bf){  // get_specific_overlap might use binary search so we have to sort
            sort_paf_buffers_by_targetID(sources, n_read);
        }
    }

    // cleanup
    free(aux.n);
    free(aux.m);
    for (int i=0; i<n_cpu; i++){
        free(aux.is_hard_add[i]);
        free(aux.is_hard_add_buf[i]);
        free(aux.buf1[i]);
        free(aux.buf2[i]);
    }
    free(aux.is_hard_add);
    free(aux.is_hard_add_buf);
    free(aux.buf1); free(aux.buf2);

    if(VERBOSE >= 0)
    {
        fprintf(stderr, "[M::%s] typeA %" PRIu64 " B %" PRIu64 ", used %.1fs\n\n", 
                        __func__, sancheckA, sancheckB, Get_T()-startTime);
        fflush(stderr);
    }
}


void get_end_match_length(ma_hit_t* edge, UC_Read* query, UC_Read* target,
uint32_t* left, uint32_t* right)
{
    (*left) = (*right) = 0;

    uint32_t query_beg, query_end, target_beg, target_end, targetLen, oLen;
    long long i;
    targetLen = Get_READ_LENGTH(R_INF, Get_tn(*edge));
    query_beg = Get_qs(*edge);
    query_end = Get_qe(*edge)-1;
    recover_UC_Read(query, &R_INF, Get_qn(*edge));
    
    if(edge->rev)
    {
        target_beg = targetLen-(Get_te(*edge)-1)-1;
        target_end = targetLen-Get_ts(*edge)-1;
        recover_UC_Read_RC(target, &R_INF, Get_tn(*edge));
    }
    else
    {
        target_beg = Get_ts(*edge);
        target_end = Get_te(*edge)-1;
        recover_UC_Read(target, &R_INF, Get_tn(*edge));
    }

    oLen = Get_qe(*edge) - Get_qs(*edge);
    /***************************left*******************************/
    for (i = 0; i < oLen; i++)
    {
        if(query->seq[query_beg+i]!=target->seq[target_beg+i]) break;
        (*left)++;
    }
    /***************************left*******************************/

    /***************************right*******************************/
    for (i = 0; i < oLen; i++)
    {
        if(query->seq[query_end-i]!=target->seq[target_end-i]) break;
        (*right)++;
    }
    /***************************right*******************************/
    
}

void normalize_ma_hit_t_single_side_aggressive(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();

    ///long long debug_prefect = 0, debug_not_bad = 0, debug_bad = 0;

    UC_Read query, target;
    init_UC_Read(&query);
    init_UC_Read(&target);

    long long i, j, index;
    uint32_t qn, tn, queryLeftLen, queryRightLen, targetLeftLen, targetRightLen;
    long long qLen_0, qLen_1, m;
    for (i = 0; i < num_sources; i++)
    {

        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index != -1)
            {
                sources[tn].buffer[index].bl = Get_qe(sources[tn].buffer[index]) 
                                                            - Get_qs(sources[tn].buffer[index]);

                if(Get_qs(sources[i].buffer[j]) == Get_ts(sources[tn].buffer[index])
                    &&
                   Get_qe(sources[i].buffer[j]) == Get_te(sources[tn].buffer[index])
                    && 
                   Get_ts(sources[i].buffer[j]) == Get_qs(sources[tn].buffer[index])
                    &&
                   Get_te(sources[i].buffer[j]) == Get_qe(sources[tn].buffer[index])
                    &&
                   sources[i].buffer[j].rev == sources[tn].buffer[index].rev)
                {
                    ///actually two edges are same, here is just to unify el, ml, no_l_indel, bl
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_prefect++;
                    continue;
                }


                get_end_match_length(&(sources[i].buffer[j]), &query, &target,
                &queryLeftLen, &queryRightLen);

                get_end_match_length(&(sources[tn].buffer[index]), &query, &target,
                &targetLeftLen, &targetRightLen);

                ///query is right
                if(queryLeftLen>0 && queryRightLen>0 && (targetLeftLen==0 || targetRightLen==0))
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_not_bad++;
                    continue;
                }


                ///target is right
                if(targetLeftLen>0 && targetRightLen>0 && (queryLeftLen==0 || queryRightLen==0))
                {
                    set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_not_bad++;
                    continue;
                }


                ///query is right
                if((queryLeftLen+queryRightLen)>(targetLeftLen+targetRightLen))
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                }///target is right
                else if((queryLeftLen+queryRightLen)<(targetLeftLen+targetRightLen))
                {
                    set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));
                }
                else
                {
                    qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                    qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);
                    if(qLen_0 >= qLen_1)
                    {
                        set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    }
                    else
                    {
                        set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));
                    }
                }

                sources[i].buffer[m] = sources[i].buffer[j];
                m++;

                ///debug_bad++;
            }
        }

        sources[i].length = m;
    }

    destory_UC_Read(&query);
    destory_UC_Read(&target);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
    }

    // fprintf(stderr, "debug_prefect: %ld, debug_not_bad: %ld, debug_bad: %ld\n", debug_prefect, debug_not_bad,
    // debug_bad);
}


void normalize_ma_hit_t_single_side(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();

    ///long long debug_prefect = 0, debug_not_bad = 0, debug_bad = 0;

    UC_Read query, target;
    init_UC_Read(&query);
    init_UC_Read(&target);

    long long i, j, index;
    uint32_t qn, tn, queryLeftLen, queryRightLen, targetLeftLen, targetRightLen;
    long long qLen_0, qLen_1, m;
    for (i = 0; i < num_sources; i++)
    {

        m = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            sources[i].buffer[j].bl = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);

            index = get_specific_overlap(&(sources[tn]), tn, qn);


            if(index != -1)
            {
                sources[tn].buffer[index].bl = Get_qe(sources[tn].buffer[index]) 
                                                            - Get_qs(sources[tn].buffer[index]);

                if(Get_qs(sources[i].buffer[j]) == Get_ts(sources[tn].buffer[index])
                    &&
                   Get_qe(sources[i].buffer[j]) == Get_te(sources[tn].buffer[index])
                    && 
                   Get_ts(sources[i].buffer[j]) == Get_qs(sources[tn].buffer[index])
                    &&
                   Get_te(sources[i].buffer[j]) == Get_qe(sources[tn].buffer[index])
                    &&
                   sources[i].buffer[j].rev == sources[tn].buffer[index].rev)
                {
                    ///actually two edges are same, here is just to unify el, ml, no_l_indel, bl
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_prefect++;
                    continue;
                }


                get_end_match_length(&(sources[i].buffer[j]), &query, &target,
                &queryLeftLen, &queryRightLen);

                get_end_match_length(&(sources[tn].buffer[index]), &query, &target,
                &targetLeftLen, &targetRightLen);

                ///query is right
                if(queryLeftLen>0 && queryRightLen>0 && (targetLeftLen==0 || targetRightLen==0))
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_not_bad++;
                    continue;
                }


                ///target is right
                if(targetLeftLen>0 && targetRightLen>0 && (queryLeftLen==0 || queryRightLen==0))
                {
                    set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));

                    sources[i].buffer[m] = sources[i].buffer[j];
                    m++;

                    ///debug_not_bad++;
                    continue;
                }

                /**
                ///query is right
                if((queryLeftLen+queryRightLen)>(targetLeftLen+targetRightLen))
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                }///target is right
                else if((queryLeftLen+queryRightLen)<(targetLeftLen+targetRightLen))
                {
                    set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));
                }
                else
                {
                    qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                    qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);
                    if(qLen_0 >= qLen_1)
                    {
                        set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                    }
                    else
                    {
                        set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));
                    }
                }
                **/



                qLen_0 = Get_qe(sources[i].buffer[j]) - Get_qs(sources[i].buffer[j]);
                qLen_1 = Get_qe(sources[tn].buffer[index]) - Get_qs(sources[tn].buffer[index]);
                if(qLen_0 >= qLen_1)
                {
                    set_reverse_overlap(&(sources[tn].buffer[index]), &(sources[i].buffer[j]));
                }
                else
                {
                    set_reverse_overlap(&(sources[i].buffer[j]), &(sources[tn].buffer[index]));
                }


                sources[i].buffer[m] = sources[i].buffer[j];
                m++;

                ///debug_bad++;
            }
        }

        sources[i].length = m;
    }

    destory_UC_Read(&query);
    destory_UC_Read(&target);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
    }

    // fprintf(stderr, "debug_prefect: %ld, debug_not_bad: %ld, debug_bad: %ld\n", debug_prefect, debug_not_bad,
    // debug_bad);
}



void drop_edges_by_trio(ma_hit_t_alloc* sources, long long num_sources)
{
    double startTime = Get_T();

    long long i, j;
    uint32_t qn, tn;
    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            if(sources[i].buffer[j].del) continue;
            
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            if(R_INF.trio_flag[qn] != AMBIGU && R_INF.trio_flag[tn] != AMBIGU)
            {
                if(R_INF.trio_flag[qn] != R_INF.trio_flag[tn])
                {
                    sources[i].buffer[j].del = 1;
                    continue;
                }
            }

            sources[i].buffer[j].del = 0;
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
    }
}

ma_hit_t* get_specific_overlap_with_del_linearsearch(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
uint32_t qn, uint32_t tn){
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i;
    for (i = 0; i < x->length; i++)
    {
        if(x->buffer[i].del) continue;
        if(coverage_cut[Get_qn(x->buffer[i])].del) continue;
        if(coverage_cut[Get_tn(x->buffer[i])].del) continue;

        if(Get_tn(x->buffer[i])==tn && Get_qn(x->buffer[i])==qn)
        {
            return &(x->buffer[i]);
        }        
    }

    return NULL;
}

ma_hit_t* get_specific_overlap_with_del_binarysearch(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
uint32_t qn, uint32_t tn){
    ma_hit_t_alloc *x = &sources[qn];
    if (x->length==0) return NULL;

    if (!x->buffer_sorted_by_tn) {
        radix_sort_hit_tn(x->buffer, x->buffer+x->length);
        x->buffer_sorted_by_tn = 1;
    }

    int low=0, high=x->length-1, idx;
    if (tn<x->buffer[0].tn || tn>x->buffer[high].tn) return NULL; 
    while (low<=high){
        idx = (high-low)/2 + low;
        if (tn<x->buffer[idx].tn){
            high = idx-1;
        }else if (tn>x->buffer[idx].tn){
            low = idx+1;
        }else{
            assert( ((uint32_t)(x->buffer[idx].qns>>32)) == qn  );
            if(x->buffer[idx].del) return NULL;
            if(coverage_cut[Get_qn(x->buffer[idx])].del) return NULL;
            if(coverage_cut[Get_tn(x->buffer[idx])].del) return NULL;
            return &(x->buffer[idx]);
        }
    }

    return NULL;
}

ma_hit_t* get_specific_overlap_with_del(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
uint32_t qn, uint32_t tn)
{
    // source but ignore entires with deleted tag; returns a handle
    if(coverage_cut[qn].del || coverage_cut[tn].del) return NULL;

    if ( (!asm_opt.get_specific_overlap_is_use_bf) || sources[qn].length<LINEAR_BF_BREAKEVEN_POINT){
		return get_specific_overlap_with_del_linearsearch(sources, coverage_cut, qn, tn);
	}else{
		return get_specific_overlap_with_del_binarysearch(sources, coverage_cut, qn, tn);
	}
}



void delete_single_edge(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, uint32_t qn, uint32_t tn)
{
    ma_hit_t* tmp = get_specific_overlap_with_del(sources, coverage_cut, qn, tn);
    if(tmp != NULL) tmp->del = 1;
}

void delete_single_edge_both_dir(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, uint32_t qn, uint32_t tn)
{
    delete_single_edge(sources, coverage_cut, qn, tn);
    delete_single_edge(sources, coverage_cut, tn, qn);
}

void delete_all_edges(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, uint32_t qn)
{
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i;
    for (i = 0; i < x->length; i++)
    {
        x->buffer[i].del = 1;
        delete_single_edge(sources, coverage_cut, Get_tn(x->buffer[i]), Get_qn(x->buffer[i]));
    }
    coverage_cut[qn].del = 1;
}


uint32_t get_real_sources_length(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp, uint32_t query)
{
    uint32_t qn = query>>1;
    if(coverage_cut[qn].del) return 0;
    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i, occ = 0;
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        if(h->del) continue;
        ///now the edge has not been removed
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        if(sq->del || st->del) continue;
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        ///if it is a contained read
        if(r < 0) continue;
        if((t.ul>>32) == query) occ++;
    }

    return occ;
}


uint32_t delete_all_edges_carefully(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp, uint32_t qn)
{
    if(coverage_cut[qn].del) return 0;

    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i, keep_edge = 0;
    uint32_t flag[2];
    flag[0] = flag[1] = 0;
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        if(h->del) continue;
        ///now the edge has not been removed
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        if(sq->del || st->del)
        {
            h->del = 1;
            delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            continue;
        }
        
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        ///if it is a contained read
        if(r < 0)
        {
            h->del = 1;
            delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            continue;
        }
        
        if(get_real_sources_length(sources, coverage_cut, max_hang, min_ovlp, (t.v^1))==1)
        {
            flag[((t.ul>>32)^1)&1] = 1;
            keep_edge++;
            continue;
        }

        // h->del = 1;
        // delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
    }


    if(keep_edge != 0)
    {
        for (i = 0; i < x->length; i++)
        {
            h = &(x->buffer[i]);
            if(h->del) continue;
            ///now the edge has not been removed
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            if(sq->del || st->del) continue;
            
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                    asm_opt.max_hang_rate, min_ovlp, &t);
            ///if it is a contained read
            if(r < 0) continue;

            if(flag[(t.ul>>32)&1] == 1) continue;
            
            if(get_real_sources_length(sources, coverage_cut, max_hang, min_ovlp, (t.v^1))==1) continue;
            
            h->del = 1;
            delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
        }
    }

    if(keep_edge == 0) coverage_cut[qn].del = 1;
    return keep_edge;
}


int ma_hit_contained_advance(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    ///uint32_t qn_num = 0, no_fully_qn_num = 0, tn_num = 0, no_fully_tn_num = 0;
    double startTime = Get_T();
    int verbose = 0;
    int ret0 = 0;
    int ret = 0;
	int32_t r;
	long long i, j, m;
	asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;

    //int dbgprint, dbgprint2;

    for (i = 0; i < n_read; ++i) 
    {
        //dbgprint = (i==38032)  || (i==67310) || (i==120968);
        //if (dbgprint) {
        //    fprintf(stderr, "[dbg::%s] saw rID qn=%d, printing all targets:\n", 
        //            __func__, (int)i);
        //    for (j=0; j<sources[i].length; j++){
        //        uint32_t tmptn = Get_tn(sources[i].buffer[j]);
        //        fprintf(stderr, "[dbg::%s] tn %d covcutdel %d strong %d el %d ovlpdel %d\n", __func__, 
        //            (int)tmptn, (int)coverage_cut[tmptn].del, 
        //            (int)sources[i].buffer[j].ml, 
        //            (int)sources[i].buffer[j].el,
        //            (int)sources[i].buffer[j].del);
        //    }
        //}
        if(coverage_cut[i].del) continue;

        for (j = 0; j < (long long)sources[i].length; j++)
        {
            h = &(sources[i].buffer[j]);
            uint32_t tn = Get_tn(*h);
            //dbgprint2 = (tn==38032)  || (tn==67310) || (tn==120968);
            //if (dbgprint||dbgprint2) fprintf(stderr, "[dbg::%s]   ckpt1 j=%d qn=%d tn=%d\n", __func__,
            //                    (int)j, (int)i, (int)Get_tn(*h));

            //check the corresponding two reads 
		    sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            /****************************may have trio bugs********************************/
            if(sq->del || st->del) continue;
            if(h->del) continue;
            /****************************may have trio bugs********************************/
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
            //if (dbgprint || dbgprint2) fprintf(stderr, "[dbg::%s]   ckpt2, r=%d\n", __func__, (int)r);
            ///r could not be MA_HT_SHORT_OVLP or MA_HT_INT
            if (r == MA_HT_QCONT) 
            {
                h->del = 1;
                delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
        
                delete_all_edges(sources, coverage_cut, Get_qn(*h));
                set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);

                ret0++;
                //if (dbgprint || dbgprint2)fprintf(stderr, "[dbg::%s] QCONT qn %d tn %d\n", __func__, 
                //        (int)Get_qn(*h), (int)Get_tn(*h));

            }
		    else if (r == MA_HT_TCONT) 
            {
                h->del = 1;
                delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));

                delete_all_edges(sources, coverage_cut, Get_tn(*h));
                set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);

                //if (dbgprint || dbgprint2) fprintf(stderr, "[dbg::%s] TCONT qn %d tn %d\n", __func__, 
                //        (int)Get_qn(*h), (int)Get_tn(*h));
            }
        }
    }

    fprintf(stderr, "[debug::%s] # QCONT: %d, used %.2f s\n", __func__, ret0, Get_T()-startTime);

    transfor_R_to_U(ruIndex);


    for (i = 0; i < n_read; ++i) 
    {
        m = 0;
        for (j = 0; j < (long long)sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            if(h->del) continue;
            ///both the qn and tn have not been deleted
            if(coverage_cut[Get_qn(*h)].del != 1 && coverage_cut[Get_tn(*h)].del != 1)
            {
                h->del = 0;
                m++;
            }
            else
            {
                h->del = 1;
            }
        }

        ///if sources[i].length == 0, that means all overlapped reads with read i are the contained reads
        if(m == 0)
        {
            fprintf(stderr, "[dbg::%s] del %d\n", __func__, (int)i);
            ret++;
            coverage_cut[i].del = 1;
        }
    }

    if(VERBOSE >= 0)
    {
        fprintf(stderr, "[M::%s] dropped %d reads, used total of %0.2f s\n\n", __func__, ret, Get_T()-startTime);
        fflush(stderr);
    }
    return ret;
}

typedef struct {
    ma_hit_t_alloc* sources;
    ma_sub_t *coverage_cut;
    uint32_t **qn;
    uint32_t **tn;
    ma_hit_t*** h_single;  // for deleting individual edges
    uint8_t **which;
    uint32_t *n;
    uint32_t *m;
    uint32_t i_batch_start;
    long long n_read;

    // A linear buffer of different length than the aboved ones, 
    //   used for deleting multiple edges.
    uint32_t *n2;
    uint32_t *m2;
    ma_hit_t*** h_multi;  // nb thread arraies of array of pointers
    uint32_t **n2_index;  // nb thread arraies of indices used to segment h_multi; size follow n, not n2
}aux_ma_hit_contain_t;

static void hamt_threaded_ma_hit_contained_advance_worker1(void *data, long i_r_, int tid){
    aux_ma_hit_contain_t *d = (aux_ma_hit_contain_t*) data;
    ma_hit_t_alloc* sources = d->sources;
    ma_sub_t *coverage_cut = d->coverage_cut;

    long long i, j, m;
    i = i_r_ + d->i_batch_start;

    if (i>=d->n_read) return;
    if (coverage_cut[i].del) return;

    ma_hit_t *h = NULL;
    ma_hit_t *tmp;
    ma_hit_t_alloc *tmp_alloc;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    asg_arc_t t;
    int32_t r;

    for (j = 0; j < (long long)sources[i].length; j++)
    {
        if (d->n[tid]==d->m[tid]){
            d->m[tid] = d->m[tid] + (d->m[tid]>>1);
            d->qn[tid] = (uint32_t*)realloc(d->qn[tid], sizeof(uint32_t) * d->m[tid]);
            d->tn[tid] = (uint32_t*)realloc(d->tn[tid], sizeof(uint32_t) * d->m[tid]);
            d->which[tid] = (uint8_t*)realloc(d->which[tid], sizeof(uint8_t) * d->m[tid]);
            d->h_single[tid] = (ma_hit_t**)realloc(d->h_single[tid], sizeof(ma_hit_t*) * d->m[tid]);
            d->n2_index[tid] = (uint32_t*)realloc(d->n2_index[tid], sizeof(uint32_t) * d->m[tid]);
        }

        h = &(sources[i].buffer[j]);
        //check the corresponding two reads 
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        /****************************may have trio bugs********************************/
        if(sq->del || st->del) continue;
        if(h->del) continue;
        /****************************may have trio bugs********************************/
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        ///r could not be MA_HT_SHORT_OVLP or MA_HT_INT
        if (r == MA_HT_QCONT) 
        {
            h->del = 1;
            // delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            // delete_all_edges(sources, coverage_cut, Get_qn(*h));
            // set_R_to_U(ruIndex, Get_qn(*h), Get_tn(*h), 0);

            // (in place of delete_single_edge)
            tmp = get_specific_overlap_with_del(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            d->h_single[tid][d->n[tid]] = tmp;
            // (in place of delete_all_edges)
            tmp_alloc = &sources[Get_qn(*h)];
            if (d->n[tid]==0){
                d->n2_index[tid][d->n[tid]] = 0;
            }else{
                d->n2_index[tid][d->n[tid]] =  d->n2_index[tid][d->n[tid]-1]+tmp_alloc->length;
            }
            if (tmp_alloc->length+d->n2[tid]>=d->m2[tid]){
                d->m2[tid] = d->m2[tid] + (d->m2[tid]>>1);
                d->h_multi[tid] = (ma_hit_t**)realloc(d->h_multi[tid], sizeof(ma_hit_t*) * d->m2[tid]);
            }
            for (int i_tmp=0; i_tmp<tmp_alloc->length; i_tmp++){
                tmp_alloc->buffer[i_tmp].del = 1;
                tmp = get_specific_overlap_with_del(sources, coverage_cut, 
                                                    Get_tn(tmp_alloc->buffer[i_tmp]), 
                                                    Get_qn(tmp_alloc->buffer[i_tmp]));
                d->h_multi[tid][d->n2[tid]] = tmp;
                d->n2[tid]++;
            }
            coverage_cut[Get_qn(*h)].del = 1;
            // (in place of set_R_to_U)
            d->qn[tid][d->n[tid]] = Get_qn(*h);
            d->tn[tid][d->n[tid]] = Get_tn(*h);
            d->which[tid][d->n[tid]] = 0;  // tell caller to get rid of qn
            d->n[tid]++;
        }
        else if (r == MA_HT_TCONT) 
        {
            h->del = 1;
            // delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            // delete_all_edges(sources, coverage_cut, Get_tn(*h));
            // set_R_to_U(ruIndex, Get_tn(*h), Get_qn(*h), 0);

            // (in place of delete_single_edge)
            tmp = get_specific_overlap_with_del(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
            d->h_single[tid][d->n[tid]] = tmp;
            // (in place of delete_all_edges)

            tmp_alloc = &sources[Get_tn(*h)];
            if (d->n[tid]==0){
                d->n2_index[tid][d->n[tid]] = 0;
            }else{
                d->n2_index[tid][d->n[tid]] =  d->n2_index[tid][d->n[tid]-1]+tmp_alloc->length;
            }
            if (tmp_alloc->length+d->n2[tid]>=d->m2[tid]){
                d->m2[tid] = d->m2[tid] + (d->m2[tid]>>1);
                d->h_multi[tid] = (ma_hit_t**)realloc(d->h_multi[tid], sizeof(ma_hit_t*) * d->m2[tid]);
            }
            for (int i_tmp=0; i_tmp<tmp_alloc->length; i_tmp++){
                tmp_alloc->buffer[i_tmp].del = 1; 
                tmp = get_specific_overlap_with_del(sources, coverage_cut, 
                                                    Get_tn(tmp_alloc->buffer[i_tmp]), 
                                                    Get_qn(tmp_alloc->buffer[i_tmp]));
                d->h_multi[tid][d->n2[tid]] = tmp;
                d->n2[tid]++;
            }
            coverage_cut[Get_tn(*h)].del = 1;
            // (in place of set_R_to_U)
            d->qn[tid][d->n[tid]] = Get_qn(*h);
            d->tn[tid][d->n[tid]] = Get_tn(*h);
            d->which[tid][d->n[tid]] = 1;  // tell caller to get rid of tn
            d->n[tid]++;
        }
    }

}
static void hamt_threaded_ma_hit_contained_advance_worker2(void *data, long i, int tid){
    aux_ma_hit_contain_t *d = (aux_ma_hit_contain_t*) data;
    ma_hit_t_alloc* sources = d->sources;
    ma_sub_t *coverage_cut = d->coverage_cut;

    long long m = 0;
    for (long long j = 0; j < (long long)sources[i].length; j++)
    {
        ma_hit_t *h = &(sources[i].buffer[j]);
        if(h->del) continue;
        ///both the qn and tn have not been deleted
        if(coverage_cut[Get_qn(*h)].del != 1 && coverage_cut[Get_tn(*h)].del != 1)
        {
            h->del = 0;
            m++;
        }
        else
        {
            h->del = 1;
        }
    }

    ///if sources[i].length == 0, that means all overlapped reads with read i are the contained reads
    if(m == 0)
    {
        d->n[tid]++;  // recycle this variable just for reporting.
        coverage_cut[i].del = 1;
    }
}

int hamt_threaded_ma_hit_contained_advance(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, 
                                           R_to_U* ruIndex){
    // FUNC
    //     Threaded ma_hit_contained_advance
    double startTime = Get_T();
    double startTime_tmp;
    double accTime_a = 0;
    double accTime_b = 0;
    int n_cpu = asm_opt.thread_num;
    int ret0 = 0; // sancheck

    // prepare aux data
    aux_ma_hit_contain_t aux;
    aux.coverage_cut = coverage_cut;
    aux.sources = sources;
    aux.n_read = n_read;
    aux.n = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.m = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.qn = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));
    aux.tn = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));
    aux.h_single = (ma_hit_t***)calloc(n_cpu, sizeof(ma_hit_t**));
    aux.n2_index = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));
    aux.which = (uint8_t**)calloc(n_cpu, sizeof(uint8_t*));  // 0 means qn of shall be removed, 1 means tn
    aux.i_batch_start = 0;
    for (int i=0; i<n_cpu; i++){
        aux.m[i] = 1024;
        aux.qn[i] = (uint32_t*)calloc(1024, sizeof(uint32_t));
        aux.tn[i] = (uint32_t*)calloc(1024, sizeof(uint32_t));
        aux.h_single[i] = (ma_hit_t**)calloc(1024, sizeof(ma_hit_t*));
        aux.n2_index[i] = (uint32_t*)calloc(1024, sizeof(uint32_t));
        aux.which[i] = (uint8_t*)calloc(1024, sizeof(uint8_t));
    }
    // (the other linear buffer)
    aux.n2 = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.m2 = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.h_multi = (ma_hit_t***)calloc(n_cpu, sizeof(ma_hit_t**));
    for (int i=0; i<n_cpu; i++){
        aux.m2[i] = 1024;
        aux.h_multi[i] = (ma_hit_t**)calloc(1024, sizeof(ma_hit_t*));
    }

    // prepare batch splits
    int batch_size = n_cpu*64;
    long long nb_batch =  n_read/((long long)batch_size)+1;
    uint32_t qn, tn;

    // step 1
    for (long long i_batch=0; i_batch<nb_batch; i_batch++){
        startTime_tmp = Get_T();

        // collect
        kt_for(n_cpu, hamt_threaded_ma_hit_contained_advance_worker1, &aux, (long)batch_size);
        accTime_a += Get_T() - startTime_tmp;

        // treat
        startTime_tmp = Get_T();
        for (int i=0; i<n_cpu; i++){
            for (int j=0; j<aux.n2[i]; j++){// in place of delete_all_edges
                if (aux.h_multi[i][j]!=NULL) {aux.h_multi[i][j]->del = 1;}
            }
            for (int j=0; j<aux.n[i]; j++){
                ret0++;
                qn = aux.qn[i][j];
                tn = aux.tn[i][j];
                if (aux.h_single[i][j]!=NULL) {aux.h_single[i][j]->del = 1;} // in place of delete_single_edge
                if (aux.which[i][j]==0){  // 
                    set_R_to_U(ruIndex, qn, tn, 0);
                }else{
                    set_R_to_U(ruIndex, tn, qn, 0);
                }
            }
        }
        accTime_b += Get_T() - startTime_tmp;

        // reset buffers
        for (int i=0; i<n_cpu; i++){
            aux.n[i] = 0;
            aux.n2[i] = 0;
        }

        // update offset
        aux.i_batch_start += batch_size;
    }

    // step 1.5: intermediate
    transfor_R_to_U(ruIndex);
    if (VERBOSE>=1) {
        fprintf(stderr, "[debug::%s] ret0 is %d\n", __func__, ret0);
        fprintf(stderr, "[prof::%s] step1 finished, part1 took %.2fs, part2 took %.2fs\n", 
                                __func__, accTime_a, accTime_b);
        
    }

    // step 2
    startTime_tmp = Get_T();
    kt_for(n_cpu, hamt_threaded_ma_hit_contained_advance_worker2, &aux, (long)n_read);
    uint64_t ret = 0;
    for (int i=0; i<n_cpu; i++){
        ret += aux.n[i];
    }

    // cleanup
    free(aux.n);
    free(aux.m);
    free(aux.n2);
    free(aux.m2);
    for (int i=0; i<n_cpu; i++){
        free(aux.qn[i]);
        free(aux.tn[i]);
        free(aux.h_single[i]);
        free(aux.n2_index[i]);
        free(aux.which[i]);
        free(aux.h_multi[i]);
    }
    free(aux.qn);
    free(aux.tn);
    free(aux.h_single);
    free(aux.n2_index);
    free(aux.h_multi);
    free(aux.which);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[prof::%s] step2 finished, took %.2fs\n", 
                                __func__, Get_T()-startTime_tmp);
        fprintf(stderr, "[M::%s] dropped %" PRIu64 " reads, takes %0.2f s\n\n", __func__, ret, Get_T()-startTime);
    }
    return ret;
}


void ma_hit_flt(ma_hit_t_alloc* sources, long long n_read, ma_sub_t *coverage_cut, int max_hang, int min_ovlp)
{
    double startTime = Get_T();
    int nb_treatedA = 0, nb_treatedB = 0;
	long long i, j, rLen;
	asg_arc_t t;

    for (i = 0; i < n_read; ++i) 
    {
        rLen = 0;
        for (j = 0; j < (long long)sources[i].length; j++)
        {
            ma_hit_t *h = &(sources[i].buffer[j]);
            if(h->del) continue;
            //check the corresponding two reads 
		    const ma_sub_t *sq = &(coverage_cut[Get_qn(*h)]);
            const ma_sub_t *st = &(coverage_cut[Get_tn(*h)]);
            int r;
		    if (sq->del || st->del) continue;
            ///[sq->s, sq->e) and [st->s, st->e) are the high coverage region in query and target
            ///here just exculde the overhang?
            ///in miniasm the 5-th option is 0.5, instead of 0.8
            /**note!!! h->qn and h->qs have been normalized by sq->s
             * h->ts and h->tn have been normalized by sq->e
             **/  
            ///here the max_hang = 1000, asm_opt.max_hang_rate = 0.8, min_ovlp = 50
            ///for me, there should not have any overhang..so r cannot be equal to MA_HT_INT
            ///sq->e - sq->s = the length of query; st->e - st->s = the length od target
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);


            ///for me, there should not have any overhang..so r cannot be equal to MA_HT_INT
            ///and I think if we use same min_ovlp in all functions, r also cannot be MA_HT_SHORT_OVLP
            ///so it does not matter we have ma_hit2arc or not
		    if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
            {
                h->del = 0;
                rLen++;
            }
            else
            {
                h->del = 1;
                delete_single_edge(sources, coverage_cut, Get_tn(*h), Get_qn(*h));
                nb_treatedA++;
            }
        
            
        }

        if(rLen == 0)
        {
            (coverage_cut)[i].del = 1;
            nb_treatedB++;
        }
    }

    fprintf(stderr, "[M::%s] typeA %d , typeB %d, used %.1fs\n", 
                    __func__, nb_treatedA, nb_treatedB, Get_T()-startTime);

}





///a is the overlap vector, n is the length of overlap vector
///min_dp is used for coverage droping
///select reads with coverage >= min_dp
void ma_hit_sub(int min_dp, ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    // min_dp is defaulted to zero (cli: `min_overlap_coverage`)
    double startTime = Get_T();

    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

	uint64_t i, j;
    int n_remained = 0, n_deleted = 0;
	kvec_t(uint32_t) b = {0,0,0};

	///all overlaps in vector a has been sorted by qns
	///so for overlaps of one reads, it must be contiguous
	for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        if(min_dp <= 1)
        {
            (*coverage_cut)[i].s = 0;
            (*coverage_cut)[i].e = readLen[i];
            (*coverage_cut)[i].del = 0;
            ++n_remained;
            continue;
        }


        kv_resize(uint32_t, b, sources[i].length);
        b.n = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            if(sources[i].buffer[j].del) continue;

            uint32_t qs, qe;
            qs = Get_qs(sources[i].buffer[j]);
            qe = Get_qe(sources[i].buffer[j]);
            kv_push(uint32_t, b, qs<<1);
			kv_push(uint32_t, b, qe<<1|1);
        }

        ///we can identify the qs and qe by the 0-th bit
		ks_introsort_uint32_t(b.n, b.a);
        ma_sub_t max, max2;
        max.s = max.e = max.del = max2.s = max2.e = max2.del = 0;
        int dp, start = 0;
        ///max is the longest subregion, max2 is the second longest subregion
        for (j = 0, dp = 0; j < b.n; ++j) 
        {
            int old_dp = dp;
            ///if a[j] is qe
            if (b.a[j]&1) 
            {
                --dp;
            }
            else
            {
                ++dp;
            } 
            
            /**
            min_dp is the coverage drop threshold
            there are two cases: 
                1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
            if one read has multiple separate sub-regions with coverage >= min_dp, 
            does miniasm only select the longest one?  // yes
            **/
            if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
            { 
                ///case 2, a[j] is qs
                start = b.a[j]>>1;
            } 
            else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
            {
                int len = (b.a[j]>>1) - start;
                if (len > (int)(max.e - max.s)) 
                {
                    max2 = max; 
                    max.s = start;
                    max.e = b.a[j]>>1;
                }
                else if (len > int(max2.e - max2.s)) 
                {
                    max2.s = start; 
                    max2.e = b.a[j]>>1;
                }
            }
        }


        ///max.e - max.s is the 
        if (max.e - max.s > 0) 
        {
            (*coverage_cut)[i].s = max.s;
            (*coverage_cut)[i].e = max.e;
            (*coverage_cut)[i].del = 0;
            ++n_remained;
		} 
        else 
        {
            (*coverage_cut)[i].s = (*coverage_cut)[i].e = 0;

            (*coverage_cut)[i].del = 1;
            n_deleted++;
        }
	}

	free(b.a);
    
    fprintf(stderr, "[M::%s] remained %d, deleted %d, used %.1fs\n", __func__, 
                n_remained, n_deleted, Get_T()-startTime);
}




int boundary_verify(uint32_t x_interval_s, uint32_t x_interval_e, ma_hit_t* map, 
char* x_buffer, char* y_buffer, All_reads* R_INF)
{
    uint32_t xs, ys, dir, x_id, y_id, x_interval_Len, y_interval_Len, y_interval_s, y_interval_e;
    dir = (*map).rev;
    xs = Get_qs((*map));
    x_id = Get_qn((*map));
    y_id = Get_tn((*map));
    long long yLen = Get_READ_LENGTH((*R_INF), y_id);

    if(dir == 1)
    {
        ys = yLen - (Get_te((*map)) - 1) - 1;         
    }
    else
    {
        ys = Get_ts((*map));
    }
    ///[x_interval_s, x_interval_e)
    x_interval_Len = x_interval_e - x_interval_s;

    ///[y_interval_s, y_interval_e]
    y_interval_s = (x_interval_s - xs) + ys;
    if(y_interval_s >= yLen)
    {
        return 0;
    }
    y_interval_e = y_interval_s + x_interval_Len - 1;
    if(y_interval_e >= yLen)
    {
        y_interval_e = yLen - 1;
    }

    if(y_interval_e < y_interval_s)
    {
        return 0;
    }
    
    y_interval_Len = y_interval_e - y_interval_s + 1;


    if(y_interval_Len <= WINDOW)
    {
        return verify_single_window(y_interval_s, y_interval_e, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF);
    }
    else
    {
        
        if(verify_single_window(y_interval_s, y_interval_s + WINDOW - 1, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF) == 0)
        {
            return 0;
        }

        if(verify_single_window(y_interval_e - WINDOW + 1, y_interval_e, ys, xs, y_id, x_id, 
                                                        dir, y_buffer, x_buffer, R_INF) == 0)
        {
            return 0;
        }

        return 1;
    }
}


void collect_sides_trio(ma_hit_t_alloc* paf, uint64_t rLen, ma_sub_t* max_left, ma_sub_t* max_right,
uint32_t trio_flag)
{
    long long j;
    uint32_t qs, qe;
    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;

        if(R_INF.trio_flag[Get_tn(paf->buffer[j])] == trio_flag) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);

        ///overlaps from left side
        if(qs == 0)
        {
            if(qs < max_left->s) max_left->s = qs;
            if(qe > max_left->e) max_left->e = qe;
        }

        ///overlaps from right side
        if(qe == rLen)
        {
            if(qs < max_right->s) max_right->s = qs;
            if(qe > max_right->e) max_right->e = qe;
        }

        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }
}


void collect_contain_trio(ma_hit_t_alloc* paf1, ma_hit_t_alloc* paf2, uint64_t rLen, 
ma_sub_t* max_left, ma_sub_t* max_right, float overlap_rate, uint32_t trio_flag)
{
    long long j, new_left_e, new_right_s;
    new_left_e = max_left->e;
    new_right_s = max_right->s;
    uint32_t qs, qe;
    ma_hit_t_alloc* paf;

    if(paf1 != NULL)
    {
        paf = paf1;
        for (j = 0; j < paf->length; j++)
        {
            if(paf->buffer[j].del) continue;

            if(R_INF.trio_flag[Get_tn(paf->buffer[j])] == trio_flag) continue;

            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }

    if(paf2 != NULL)
    {
        paf = paf2;
        for (j = 0; j < paf->length; j++)
        {
            if(paf->buffer[j].del) continue;

            if(R_INF.trio_flag[Get_tn(paf->buffer[j])] == trio_flag) continue;
            
            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }
    

    max_left->e = new_left_e;
    max_right->s = new_right_s;
}


void collect_sides(ma_hit_t_alloc* paf, uint64_t rLen, ma_sub_t* max_left, ma_sub_t* max_right)
{
    long long j;
    uint32_t qs, qe;
    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);


        ///overlaps from left side
        if(qs == 0)
        {
            if(qs < max_left->s) max_left->s = qs;
            if(qe > max_left->e) max_left->e = qe;
        }

        ///overlaps from right side
        if(qe == rLen)
        {
            if(qs < max_right->s) max_right->s = qs;
            if(qe > max_right->e) max_right->e = qe;
        }

        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }
}

uint64_t collect_sides_also_count(ma_hit_t_alloc* paf, uint64_t rLen, ma_sub_t* max_left, ma_sub_t* max_right)
{
    long long j;
    uint32_t qs, qe;
    uint64_t bp_left=0, bp_right=0, cnt_left=0, cnt_right=0;
    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);


        ///overlaps from left side
        if(qs == 0)
        {
            if(qs < max_left->s) max_left->s = qs;
            if(qe > max_left->e) max_left->e = qe;
            bp_left+=qe-qs;
        }

        ///overlaps from right side
        if(qe == rLen)
        {
            if(qs < max_right->s) max_right->s = qs;
            if(qe > max_right->e) max_right->e = qe;
            bp_right+=qe-qs;
        }
        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }
    
    return bp_left<<32 | bp_right;
}

void collect_contain(ma_hit_t_alloc* paf1, ma_hit_t_alloc* paf2, uint64_t rLen, 
ma_sub_t* max_left, ma_sub_t* max_right, float overlap_rate)
{
    long long j, new_left_e, new_right_s;
    new_left_e = max_left->e;
    new_right_s = max_right->s;
    uint32_t qs, qe;
    ma_hit_t_alloc* paf;

    if(paf1 != NULL)
    {
        paf = paf1;
        for (j = 0; j < paf->length; j++)
        {
            if(paf->buffer[j].del) continue;

            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }

    if(paf2 != NULL)
    {
        paf = paf2;
        for (j = 0; j < paf->length; j++)
        {
            if(paf->buffer[j].del) continue;

            qs = Get_qs(paf->buffer[j]);
            qe = Get_qe(paf->buffer[j]);
            ///check contained overlaps
            if(qs != 0 && qe != rLen)
            {
                ///[qs, qe), [max_left.s, max_left.e)
                if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
                {
                    ///if(qe > max_left->e) max_left->e = qe;
                    if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
                }

                ///[qs, qe), [max_right.s, max_right.e)
                if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
                {
                    ///if(qs < max_right->s) max_right->s = qs;
                    if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
                }
            }
        }
    }
    

    max_left->e = new_left_e;
    max_right->s = new_right_s;
}


uint32_t collect_contain_also_count(ma_hit_t_alloc* paf1, uint64_t rLen, 
                                ma_sub_t* max_left, ma_sub_t* max_right, float overlap_rate)
{
    long long j, new_left_e, new_right_s;
    new_left_e = max_left->e;
    new_right_s = max_right->s;
    uint32_t qs, qe;
    ma_hit_t_alloc* paf;

    uint32_t bp_mid = 0, cnt=0;

    paf = paf1;
    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;
        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///check contained overlaps
        if(qs != 0 && qe != rLen)
        {
            bp_mid += qe-qs;
            ///[qs, qe), [max_left.s, max_left.e)
            if(qs < max_left->e && qe > max_left->e && max_left->e - qs > (overlap_rate * (qe -qs)))
            {
                ///if(qe > max_left->e) max_left->e = qe;
                if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
            }

            ///[qs, qe), [max_right.s, max_right.e)
            if(qs < max_right->s && qe > max_right->s && qe - max_right->s > (overlap_rate * (qe -qs)))
            {
                ///if(qs < max_right->s) max_right->s = qs;
                if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
            }
        }
    }
    return bp_mid/cnt;
}


int intersection_check(ma_hit_t_alloc* paf, uint64_t rLen, uint32_t interval_s, uint32_t interval_e)
{
    long long j, cov = 0;
    uint32_t qs, qe;

    for (j = 0; j < paf->length; j++)
    {
        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///[interval_s, interval_e) must be at least contained at one of the [qs, qe)
        if(qs<=interval_s && qe>=interval_e)
        {
            cov++;
        }
    }

    return cov;
}


int intersection_check_by_base(ma_hit_t_alloc* paf, uint64_t rLen, uint32_t interval_s, uint32_t interval_e,
char* bq, char* bt)
{
    long long j;
    uint32_t qs, qe;

    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///[interval_s, interval_e) must be at least contained at one of the [qs, qe)
        if(qs<=interval_s && qe>=interval_e)
        {
            if(boundary_verify(interval_s, interval_e, &(paf->buffer[j]), bq, bt, &R_INF) == 0)
            {
                return 1;
            }
        }
    }

    return 0;
}

void print_overlaps(ma_hit_t_alloc* paf, long long rLen, long long interval_s, long long interval_e)
{
    long long j;

    fprintf(stderr, "left: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qs(paf->buffer[j]) == 0)
        {
            fprintf(stderr, "shouldn't happen; interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);

            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }   
    }

    fprintf(stderr, "right: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qe(paf->buffer[j]) == rLen)
        {
            fprintf(stderr, "shouldn't happen interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);
            
            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }
    }


    fprintf(stderr, "middle: \n");
    for (j = 0; j < paf->length; j++)
    {
        if(Get_qs(paf->buffer[j]) != 0 && Get_qe(paf->buffer[j]) != rLen)
        {
            fprintf(stderr, "shouldn't happen interval_s: %lld, interval_e: %lld, qn: %u, tn: %u, j: %lld, qs: %u, qe: %u, ts: %u, te: %u, dir: %u\n",
            interval_s, interval_e,
            Get_qn(paf->buffer[j]), Get_tn(paf->buffer[j]),
            j, Get_qs(paf->buffer[j]), Get_qe(paf->buffer[j]),
            Get_ts(paf->buffer[j]), Get_te(paf->buffer[j]),
            paf->buffer[j].rev);
            
            fprintf(stderr, "%.*s\n", (int)Get_NAME_LENGTH(R_INF, Get_tn(paf->buffer[j])),
                Get_NAME(R_INF, Get_tn(paf->buffer[j])));
        }
    }

}


void detect_chimeric_reads(ma_hit_t_alloc* paf, long long n_read, uint64_t* readLen, 
ma_sub_t* coverage_cut, float shift_rate)
{
    double startTime = Get_T();
    init_aux_table();
    long long i, rLen, /**cov,**/ n_simple_remove = 0, n_complex_remove = 0, n_complex_remove_real = 0;
    uint32_t interval_s, interval_e;
    ma_sub_t max_left, max_right;
    kvec_t(char) b_q = {0,0,0};
    kvec_t(char) b_t = {0,0,0};
    for (i = 0; i < n_read; ++i) 
    {
        coverage_cut[i].c = PRIMARY_LABLE;
        rLen = readLen[i];


        max_left.s = max_right.s = rLen;
        max_left.e = max_right.e = 0;

        collect_sides(&(paf[i]), rLen, &max_left, &max_right);
        ///collect_sides(&(rev_paf[i]), rLen, &max_left, &max_right);
        ///that means this read is an end node
        if(max_left.s == rLen || max_right.s == rLen)
        {
            continue;
        }

        collect_contain(&(paf[i]), NULL, rLen, &max_left, &max_right, 0.1);
        ///collect_contain(&(paf[i]), &(rev_paf[i]), rLen, &max_left, &max_right, 0.1);

        ////shift_rate should be (asm_opt.max_ov_diff_final*2)
        ///this read is a normal read
        if (max_left.e > max_right.s && (max_left.e - max_right.s >= rLen * shift_rate))
        {
            continue;
        }

        ///simple chimeric reads
        if(max_left.e <= max_right.s)
        {
            delete_all_edges(paf, coverage_cut, i);            
            n_simple_remove++;
            continue;
        }

        ///now max_left.e >  max_right.s && max_left.e - max_right.s is small enough
        //[interval_s, interval_e)
        interval_s = max_right.s;
        interval_e = max_left.e;

        /**
        cov = 0;
        cov += intersection_check(&(paf[i]), rLen, interval_s, interval_e);
        cov += intersection_check(&(rev_paf[i]), rLen, interval_s, interval_e);
        if(interval_e - interval_s < WINDOW && cov <= 2)
        {
            coverage_cut[i].del = 1;
            paf[i].length = 0;
            n_complex_remove_real++;
        }
        else**/
        {
            kv_resize(char, b_q, WINDOW*4+20);
            kv_resize(char, b_t, WINDOW*4+20);
            if(intersection_check_by_base(&(paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)
              /**||
              intersection_check_by_base(&(rev_paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)**/)
            {
                delete_all_edges(paf, coverage_cut, i); 
                n_complex_remove_real++;
            }
        }

        n_complex_remove++;
    }

    free(b_q.a);
    free(b_t.a);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s, n_simple_remove: %lld, n_complex_remove: %lld/%lld\n\n", 
        __func__, Get_T()-startTime, n_simple_remove, n_complex_remove_real, n_complex_remove);
    }
}


int hamt_ovlp_is_end_read(ma_hit_t_alloc* sources, ma_hit_t_alloc *reverse_sources, long long i_read, uint64_t *readLen, int overhang_allowance){
    // RETURN
    //    1 if yes
    //    0 if perfectly no
    //    -1 if imperfectly no
    //    -2 if max overhang is larger than overhang_allowance
    // NOTE
    //    Collect_sides_also_count touches only ideal hits (i.e. start is 0 or end is readlength),
    //    but here we might have to deal with erroneours reads that's not actually an end read & doesn't overlap well however.
    ma_hit_t_alloc *paf_of_read = &sources[i_read];
    ma_hit_t_alloc *paf_of_read_rev;
    if (reverse_sources!=NULL) {paf_of_read_rev = &reverse_sources[i_read];}
    uint64_t readlength = readLen[i_read]; 

    ma_sub_t max_left, max_right;
    ma_hit_t *handle;
    max_left.s = max_right.s = readlength;
    max_left.e = max_right.e = 0;
    collect_sides_also_count(paf_of_read, readlength, &max_left, &max_right);
    
    if( max_left.s==0 && max_right.e==readlength){return 0;}
    else {
        uint32_t left_passed=0, right_passed=0, t_lefthang, t_righthang;
        uint32_t qs, qe;
        uint64_t t_readlength;
        // intra-haplotype
        for (int i=0; i<paf_of_read->length; i++){
            if (paf_of_read->buffer[i].del) {continue;}
            handle = &paf_of_read->buffer[i];
            qs = (uint32_t)handle->qns;
            qe = handle->qe;
            t_readlength = readLen[handle->tn];

            if (handle->rev==0){  // "left" assumes laying the query left to right and use the orientation observed
                t_lefthang = handle->ts;
                t_righthang = t_readlength - handle->te;
            }else{
                t_lefthang = t_readlength - handle->te;
                t_righthang = handle->ts;
            }
            // fprintf(stderr, "    check target %.*s, lefthang:%d, righthang:%d\n", (int)Get_NAME_LENGTH(R_INF, handle->tn), Get_NAME(R_INF, handle->tn),
                            // t_lefthang, t_righthang);

            // sancheck
            if (qs>overhang_allowance && (readlength-qe)>overhang_allowance &&
                t_lefthang>overhang_allowance && t_righthang>overhang_allowance){
                return -2;
            }
            // check coordinates
            if (t_lefthang>=qs){left_passed = 1;}
            if (t_righthang>=(readlength-qe)) {right_passed = 1;}
        }
        // inter-haplotype
        if (reverse_sources!=NULL){
            for (int i=0; i<paf_of_read_rev->length; i++){
                if (paf_of_read_rev->buffer[i].del) {continue;}
                handle = &paf_of_read_rev->buffer[i];
                qs = (uint32_t)handle->qns;
                qe = handle->qe;
                t_readlength = readLen[handle->tn];

                if (handle->rev==0){  // "left" assumes laying the query left to right and use the orientation observed
                    t_lefthang = handle->ts;
                    t_righthang = t_readlength - handle->te;
                }else{
                    t_lefthang = t_readlength - handle->te;
                    t_righthang = handle->ts;
                }
                // fprintf(stderr, "    check target %.*s, lefthang:%d, righthang:%d\n", (int)Get_NAME_LENGTH(R_INF, handle->tn), Get_NAME(R_INF, handle->tn),
                                // t_lefthang, t_righthang);

                // sancheck
                if (qs>overhang_allowance && (readlength-qe)>overhang_allowance &&
                    t_lefthang>overhang_allowance && t_righthang>overhang_allowance){
                    return -2;
                }
                // check coordinates
                if (t_lefthang>=qs){left_passed = 1;}
                if (t_righthang>=(readlength-qe)) {right_passed = 1;}
            }
        }

        if (left_passed && right_passed){
            // this is effectively: query read maps with some other reads in the middle,
            //                      but towards the end they don't match.
            return -1;
        }
    }
    return 1;
}

int hamt_ovlp_is_contained_read(ma_hit_t_alloc* paf_of_read, long long readlength){
    ma_sub_t max_left, max_right;
    max_left.s = max_right.s = readlength;
    max_left.e = max_right.e = 0;
    collect_sides_also_count(paf_of_read, readlength, &max_left, &max_right);
    if(max_left.s == 0 && max_right.s == 0 && max_left.e==readlength && max_right.e==readlength) {return 1;}
    else {return 0;}
}
int hamt_ovlp_is_contained_read2(ma_hit_t_alloc* sources, long long i_read, uint64_t *readLen){
    return hamt_ovlp_is_contained_read(&sources[i_read], readLen[i_read]);
}

int hamt_ovlp_contains_howmany_reads(ma_hit_t_alloc* sources, long long i_read, uint64_t *readLen){
    ma_hit_t *h;
    int ret = 0;
    for (int i=0; i<sources[i_read].length; i++){
        h = &sources[i_read].buffer[i];
        if (h->del) {continue;}
        if ( (h->ts==0 && h->te==readLen[h->tn]) || (h->ts==readLen[h->tn] && h->te==0) ){
            ret+=1;
        }
    }
    return ret;
}

inline uint16_t hamt_count_error(long long i_read){
    // NOTE
    //    cigar and round2.cigar are temporary buffer during ovec, fully corrected reads will have them as all zero
    return R_INF.nb_error_corrected[i_read];
}

void detect_chimeric_reads_conservative(ma_hit_t_alloc* paf, long long n_read, uint64_t* readLen, 
ma_sub_t* coverage_cut, float shift_rate)
{
    // hamt NOTE
    //     Added conservative heuristic for the simple case, because we observed read chain breaking at
    //       spots that only supported by one read, with or without other complications,
    //       in species/strain with low coverage. 
    int verbose =0;

    double startTime = Get_T();
    init_aux_table();
    long long i, rLen, /**cov,**/ n_simple_remove = 0, n_complex_remove = 0, n_complex_remove_real = 0;
    uint32_t interval_s, interval_e;
    ma_sub_t max_left, max_right;
    kvec_t(char) b_q = {0,0,0};
    kvec_t(char) b_t = {0,0,0};

    uint64_t est_cov;  // "estimated coverage"
    uint32_t est_cov_left, est_cov_right, est_cov_mid;

    for (i = 0; i < n_read; ++i) 
    {
        coverage_cut[i].c = PRIMARY_LABLE;
        rLen = readLen[i];

        max_left.s = max_right.s = rLen;
        max_left.e = max_right.e = 0;

        est_cov = collect_sides_also_count(&(paf[i]), rLen, &max_left, &max_right);
        
        if(max_left.s == rLen || max_right.s == rLen)  // the read is a end node
        {
            continue;
        }

        collect_contain(&(paf[i]), NULL, rLen, &max_left, &max_right, 0.1);
        ///collect_contain(&(paf[i]), &(rev_paf[i]), rLen, &max_left, &max_right, 0.1);

        ////shift_rate should be (asm_opt.max_ov_diff_final*2)
        ///this read is a normal read
        if (max_left.e > max_right.s && (max_left.e - max_right.s >= rLen * shift_rate))
        {
            continue;
        }

        ///simple chimeric reads
        if(max_left.e <= max_right.s)
        {
            //////////// extra check: 
            // don't drop if 1) both sides are low coverage
            //               2) the read in question had little error before ovec
            // check estimated coverage
            est_cov_left = ((uint32_t)(est_cov>>32)) / (max_left.e-max_left.s);
            est_cov_right = (uint32_t)est_cov / (max_right.e-max_right.s);
            if (est_cov_left<=5 && est_cov_right<=5){
                if (verbose) {fprintf(stderr, "chimera: sparing read %.*s\n", (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));}
                continue;
            }

            // all clear, drop it
            delete_all_edges(paf, coverage_cut, i);            
            n_simple_remove++;
            continue;
        }

        ///now max_left.e >  max_right.s && max_left.e - max_right.s is small enough
        //[interval_s, interval_e)
        interval_s = max_right.s;
        interval_e = max_left.e;

        /**
        cov = 0;
        cov += intersection_check(&(paf[i]), rLen, interval_s, interval_e);
        cov += intersection_check(&(rev_paf[i]), rLen, interval_s, interval_e);
        if(interval_e - interval_s < WINDOW && cov <= 2)
        {
            coverage_cut[i].del = 1;
            paf[i].length = 0;
            n_complex_remove_real++;
        }
        else**/
        {
            kv_resize(char, b_q, WINDOW*4+20);
            kv_resize(char, b_t, WINDOW*4+20);
            if(intersection_check_by_base(&(paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)
              /**||
              intersection_check_by_base(&(rev_paf[i]), rLen, interval_s, interval_e, b_q.a, b_t.a)**/)
            {
                delete_all_edges(paf, coverage_cut, i); 
                n_complex_remove_real++;
            }
        }

        n_complex_remove++;
    }

    free(b_q.a);
    free(b_t.a);

    fprintf(stderr, "[M::%s] n_simple_remove: %lld, n_complex_remove: %lld/%lld, used %.1f\n", 
        __func__, n_simple_remove, n_complex_remove_real, n_complex_remove, Get_T()-startTime);
}

int hamt_hit_contained_nested_containment(ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, 
                                          ma_sub_t *coverage_cut,
                                          uint64_t *readLen){
    // FUNC
    //     Called prior to hamt_hit_contained_drop_singleton_worker.
    //     When a read is nestedly contained, check up on the parents and ditch the read (shrotest in the nest)
    //      if it has a parent of the same hap && 
    //      the parent's parents are also of the same hap as the read's parents.
    // NOTE
    //     Serial. 
    // RETURN
    //     number of reads treated.
    int ret = 0;
    int verbose = 0;

    // read's parents
    uint32_t *buf = (uint32_t*)malloc(sizeof(uint32_t)*16);
    int buf_size = 16;
    int buf_len = 0;
    // parents of parents
    uint32_t *buf_pp = (uint32_t*)malloc(sizeof(uint32_t)*16);
    int buf_pp_size = 16;
    int buf_pp_len = 0;
    // for truncating buf
    int tmp_len = 0;


    int cnt_parent = 0;
    int cnt_parent_rev = 0;

    for (uint32_t i_r=0; i_r<R_INF.total_reads; i_r++){
        if (verbose>1) fprintf(stderr, "[debug::%s] at read %.*s (%d)\n", __func__,
                        (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r), (int)i_r);
        // reset
        buf_len = 0;
        buf_pp_len = 0;
        tmp_len = 0;

        // check if read is contained
        cnt_parent = hamt_read_is_contained_by_others(sources, i_r, readLen[i_r],
                                                      &buf, &buf_size, &buf_len);
        cnt_parent_rev = hamt_read_is_contained_by_others(reverse_sources, i_r, readLen[i_r],
                                                      &buf, &buf_size, &buf_len);
        cnt_parent_rev = 0;
        if (cnt_parent==0) continue;
        if (cnt_parent_rev>0) continue;  // the read has parents from different haps, don't touch it here.
        // NOTE/potential bug:
        //     `ma_hit_contained_advance` only checks on sources, not reverse_sources
        //      there's 'containment' in the later collection. 
        //      What to do about them?

        // check if any of the parents are contained
        if (verbose>1) {fprintf(stderr, "check parent containment, buf len is %d\n", buf_len);}
        cnt_parent = 0;
        int tmp_cnt_parent;
        for (int i=0; i<buf_len; i++){
            tmp_cnt_parent = hamt_read_is_contained_by_others(sources, buf[i], readLen[buf[i]], 
                                                            &buf_pp, &buf_pp_size, &buf_pp_len);
            if (tmp_cnt_parent>0){
                buf[tmp_len] = buf[i];
                tmp_len++;
            }
            cnt_parent += tmp_cnt_parent;
        }
        if (cnt_parent==0) continue;  // not nested containment
        buf_len = tmp_len;  // now buf only contains parent reads that are also contained by some other reads.

        // check parent reads' parents, pairwise
        if (verbose>1) {fprintf(stderr, "check parent of parents\n");}
        uint8_t stat=0;
        for (int i=0; i<buf_pp_len-1; i++){
            for (int j=i+1; j<buf_pp_len; j++){
                stat = hamt_check_if_reads_overlap(buf_pp[i], buf_pp[j], sources, coverage_cut, 0);
                if (stat!=3){break;}
            }
            if (stat!=3){break;}
        }
        if (stat!=3) {continue;}

        // if we reach here, 
        // the current read can be safely removed
        if (verbose>1) {fprintf(stderr, "drop!\n");}
        coverage_cut[i_r].del = 1;
        for (int i=0; i<buf_len; i++){
            delete_single_edge_both_dir(sources, coverage_cut, buf[i], (uint32_t)i_r); 
        }
        ret+=1;
        if (verbose) {fprintf(stderr, "[debug::%s] removed %.*s (%d)\n", __func__,
                                    (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                                    (int)i_r);
                     }
    }

    free(buf);
    free(buf_pp);
    if (VERBOSE){
        fprintf(stderr, "[M::%s] removed %d reads\n", __func__, ret);
    }

    return ret;
}


typedef struct {
    ma_hit_t_alloc *sources;
    ma_hit_t_alloc *reverse_sources;
    long long n_read;
    uint64_t *readLen;
    ma_sub_t *coverage_cut;
    int counter;  // NOTE this might RACE during update, but it's an output and not used elsewhere
}hitcontain_aux_t;

static void hamt_hit_contained_drop_singleton_worker(void *data, long i_r, int tid){  // callback for kt_for()
    // NOTE / might have bugs
    //     It might be questionable to see a read to be contained by one and only one another read.
    //     One case is that the other read comes from another haplotype (phasing variant is 
    //        outside of the exact match region), and that haplotype breaks while
    //        the haplotype of this contained read has a better continuity.
    //        We would want to favor the contained reads if so.
    //     The heuristic is to disallow read containment if the shorter read is not
    //        contained by at least two haplotypes. 
    //     By the way, checking whether two reads are from the same haplotype is trivial,
    //        but getting the exact position of the phasing variants is convoluted and
    //        must be done during the ovec stage (data discarded immediately, see maybe: `worker_ovec`).
    // FUNC
    //     Drop arcs to protect certain shorter reads from being considered as contained.
    //     This function won't touch the reads. Only arcs.

    int verbose = 0;  // debug

    fprintf(stderr, "[debug::%s] > at read %.*s (%d)\n", __func__,
                    (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                     (int)i_r);
    
    // input
    hitcontain_aux_t *d = (hitcontain_aux_t*)data;
    ma_hit_t_alloc *sources = d->sources;
    ma_hit_t_alloc *reverse_sources = d->reverse_sources;
    long long n_read = d->n_read;
    uint64_t *readLen = d->readLen;
    ma_sub_t *coverage_cut = d->coverage_cut;

    // early termination
    if (coverage_cut[i_r].del) {return;}
    if (R_INF.mask_readnorm[i_r] & 1){return;}
    
    // aux
    ma_hit_t_alloc *h = &sources[i_r];
    ma_hit_t_alloc *h_rev = &reverse_sources[i_r];
    ma_hit_t *hit = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    int32_t r;  // for ma_hit2arc return value
    asg_arc_t t;  // placeholder
    int need_to_protect = 0;
    uint32_t v, w;
    int i, j;

    // buffer
    int buf_len=0;
    int buf_size = 16;
    uint32_t *buf = (uint32_t*)malloc(sizeof(uint32_t) * buf_size);
    ma_hit_t **hits = (ma_hit_t**)malloc(sizeof(ma_hit_t*) * buf_size);  // store the handles, in case we want to drop arcs


    // check containments
    for (i=0; i<h->length; i++){
        hit = &(h->buffer[i]);
        sq = &(coverage_cut[Get_qn(*hit)]);
        st = &(coverage_cut[Get_tn(*hit)]);
        if(sq->del || st->del) continue;
        if(hit->del) continue;
        r = ma_hit2arc(hit, sq->e - sq->s, st->e - st->s, 
                       asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        if (r==MA_HT_QCONT){
            buf[buf_len] = hit->tn;
            hits[buf_len] = hit;
            buf_len++;
            if (buf_len==buf_size){
                buf_size = buf_size + (buf_size>>1);
                buf = (uint32_t*)realloc(buf, sizeof(uint32_t)*buf_size);
                hits = (ma_hit_t**)realloc(hits, sizeof(ma_hit_t*) * buf_size);
            }
        }
    }
    for (i=0; i<h_rev->length; i++){
        hit = &(h_rev->buffer[i]);
        sq = &(coverage_cut[Get_qn(*hit)]);
        st = &(coverage_cut[Get_tn(*hit)]);
        if(sq->del || st->del) continue;
        if(hit->del) continue;
        r = ma_hit2arc(hit, sq->e - sq->s, st->e - st->s, 
                       asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        if (r==MA_HT_QCONT){
            buf[buf_len] = hit->tn;
            hits[buf_len] = hit;
            buf_len++;
            if (buf_len==buf_size){
                buf_size = buf_size + (buf_size>>1);
                buf = (uint32_t*)realloc(buf, sizeof(uint32_t)*buf_size);
                hits = (ma_hit_t**)realloc(hits, sizeof(ma_hit_t*) * buf_size);
            }
        }
    }

    // check up on the "parent" reads
    if (buf_len==0){  // read not contained
        need_to_protect = -1;
        goto finish;
    }else{  // check the "parent" reads pairwise-ly
        need_to_protect = 1;
        for (i=0; i<buf_len-1; i++){
            for (j=i+1; j<buf_len; j++){
                v = buf[i];
                w = buf[j];
                if (hamt_check_if_reads_overlap(v, w, reverse_sources, coverage_cut, 0)>0){
                    need_to_protect = 0;
                    goto finish;
                }
            }
        }
    }
    
finish:
    if (need_to_protect>0){  // drop arcs between the shorter read and its "parent" reads
        if (verbose){
            fprintf(stderr, "[debug::%s]     protect %.*s \n", __func__, 
                            (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r));
            fprintf(stderr, "[debug::%s]       traceback - the parent reads:\n", __func__);
            for (i=0; i<buf_len; i++){
                fprintf(stderr, "[debug::%s]         %.*s\n", __func__,
                            (int)Get_NAME_LENGTH(R_INF, hits[i]->tn), Get_NAME(R_INF, hits[i]->tn));
            }
        }
        d->counter++;
        for (i=0; i<buf_len; i++){
            delete_single_edge_both_dir(sources, coverage_cut, (uint32_t)i_r, hits[i]->tn);
            delete_single_edge_both_dir(reverse_sources, coverage_cut, (uint32_t)i_r, hits[i]->tn);
            
        }
    }else if (need_to_protect==0){
        if (verbose){
            fprintf(stderr, "[debug::%s]     contained but not protecting %.*s \n", __func__, 
                            (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r));
            fprintf(stderr, "[debug::%s]       traceback - the parent reads:\n", __func__);
            for (i=0; i<buf_len; i++){
                fprintf(stderr, "[debug::%s]         %.*s\n", __func__,
                            (int)Get_NAME_LENGTH(R_INF, hits[i]->tn), Get_NAME(R_INF, hits[i]->tn));
            }
        }
    }else{  // read not contained
        ;
    }

    free(buf);
    free(hits);
}

static void hamt_hit_contained_drop_singleton_worker_v2(void *data, long i_r, int tid){  // callback for kt_for()
    // NOTE / might have bugs
    //     v1 doesn't work if there's unfortunately many haps. Consider:
    //        (x and o stand for phasing variants)
    //        read A  ------x------------------x------            (haplotypeA)
    //        read B     ---x------------------o----------        (haplotypeB)
    //        read C             --------------o---               (    ?     )
    //        read D  ------o----------------                     (haplotypeD)
    //     In this case we have 3 haplotypes: A, B, D. The read C could belong to hapB or hapD.
    //     If the coverage of haplotype D is sufficient, it doesn't matter if 
    //       we collapsed C onto B based on containment,
    //       or based on v1's criteria (C can be removed because it has parents from at least 2 haps) .
    //     The read D will find another read to go through this region.
    //     However, if haplotype D happens to have to rely on C->D in the region,
    //        losing read C is losing haplotype D's continuity.
    //     Note that this doesn't need the region to be very low coverage overally,
    //        and the above plot is a real case - there's only 2 phasing variants.
    // NOTE
    //     If coverage is apparently sufficient, this can be harmful.
    // FUNC
    //     Drop arcs to protect certain shorter reads from being considered as contained.
    //     This function won't touch the reads. Only arcs.

    int verbose = 0;  // debug

    if (verbose>1) {fprintf(stderr, "[debug::%s] > at read %.*s (%d)\n", __func__,
                    (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                     (int)i_r);
    }
    
    // input
    hitcontain_aux_t *d = (hitcontain_aux_t*)data;
    ma_hit_t_alloc *sources = d->sources;
    ma_hit_t_alloc *reverse_sources = d->reverse_sources;
    long long n_read = d->n_read;
    uint64_t *readLen = d->readLen;
    ma_sub_t *coverage_cut = d->coverage_cut;

    // early termination
    if (coverage_cut[i_r].del) {return;}
    if (R_INF.mask_readnorm[i_r] & 1){return;}
    if (sources[i_r].length>10) {return;}
    
    // aux
    ma_hit_t_alloc *h = &sources[i_r];
    ma_hit_t_alloc *h_rev = &reverse_sources[i_r];
    ma_hit_t *hit = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    int32_t r;  // for ma_hit2arc return value
    asg_arc_t t;  // placeholder
    int need_to_protect = 0;
    uint32_t v, w;
    int i, j;
    int homo_strong_cnt = 0;

    // debug aux
    int nb_total_ovlp = 0;
    int nb_weak_ovlp = 0;  // for debug; how many weak/total non-containment reads overlap with the query read

    // buffer
    int buf_len=0;
    int buf_size = 16;
    uint32_t *buf = (uint32_t*)malloc(sizeof(uint32_t) * buf_size);
    ma_hit_t **hits = (ma_hit_t**)malloc(sizeof(ma_hit_t*) * buf_size);  // store the handles, in case we want to drop arcs


    // check containments
    // Note that we need to check both paf and reverse_paf,
    //   since read's neighbours will be compared against all possible parent haplotypes.
    // This is different from the actual containment collapsing routine where
    //   only intra-haplotype arcs are considered.
    // (intra)
    for (i=0; i<h->length; i++){
        hit = &(h->buffer[i]);
        sq = &(coverage_cut[Get_qn(*hit)]);
        st = &(coverage_cut[Get_tn(*hit)]);
        if(sq->del || st->del) continue;
        if(hit->del) continue;
        r = ma_hit2arc(hit, sq->e - sq->s, st->e - st->s, 
                       asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        if (r==MA_HT_QCONT){
            buf[buf_len] = hit->tn;
            hits[buf_len] = hit;
            buf_len++;
            if (buf_len==buf_size){
                buf_size = buf_size + (buf_size>>1);
                buf = (uint32_t*)realloc(buf, sizeof(uint32_t)*buf_size);
                hits = (ma_hit_t**)realloc(hits, sizeof(ma_hit_t*) * buf_size);
            }
        }
    }
    if (buf_len==0){goto finish;}  // read doesn't have a chance to be collapsed due to containment anyway.
    // (inter)
    for (i=0; i<h_rev->length; i++){
        hit = &(h_rev->buffer[i]);
        sq = &(coverage_cut[Get_qn(*hit)]);
        st = &(coverage_cut[Get_tn(*hit)]);
        if(sq->del || st->del) continue;
        if(hit->del) continue;
        r = ma_hit2arc(hit, sq->e - sq->s, st->e - st->s, 
                       asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        if (r==MA_HT_QCONT){
            buf[buf_len] = hit->tn;
            hits[buf_len] = hit;
            buf_len++;
            if (buf_len==buf_size){
                buf_size = buf_size + (buf_size>>1);
                buf = (uint32_t*)realloc(buf, sizeof(uint32_t)*buf_size);
                hits = (ma_hit_t**)realloc(hits, sizeof(ma_hit_t*) * buf_size);
            }
        }
    }

    // TODO
    // We don't want to protect too many contained reads, the graph complexity explodes.
    // Would be ideal if we can esitmate coverage for each haplotype, but that needs more effort.
    // The compromise here is to assume that the shorter read is safe to drop if
    //    it has a bunch of parent reads or arcs, which implies this region is 
    //    getting at least some coverage.
    if (buf_len>3){goto finish;}


    // check neighbours
    v = (uint32_t)i_r;
    for (i=0; i<h->length; i++){ 
        hit = &(h->buffer[i]);
        w = hit->tn;
        if (verbose>1){
                fprintf(stderr, "[debug::%s]     = will check neighbor %.*s\n",
                            __func__, (int)Get_NAME_LENGTH(R_INF, w), Get_NAME(R_INF, w));
            }

        sq = &(coverage_cut[Get_qn(*hit)]);
        st = &(coverage_cut[Get_tn(*hit)]);
        // if(sq->del || st->del) continue;
        // if (hit->del) continue;
        r = ma_hit2arc(hit, sq->e - sq->s, st->e - st->s, 
                       asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        // if (r<0) continue;  // ignore if the ovlp is containment(either way), short, or there's large overhang
        if (r==MA_HT_QCONT || r==MA_HT_TCONT) continue;
        nb_total_ovlp++;
        if (hit->ml){continue;}  // ignore strong overlap
        nb_weak_ovlp++;

        if (verbose>1) {fprintf(stderr, "[debug::%s]      (checkpoint)\n", __func__);}

        // check this neighbour against all parent reads
        if (hamt_read_het_to_all_given_reads(sources, w, buf, buf_len, 1)){
            // additionally, require this read to have no strong intra-haplotype arc
            homo_strong_cnt = 0;
            for (j=0; j<sources[w].length; j++){
                if (sources[w].buffer[j].ml) homo_strong_cnt++;
            }
            if (homo_strong_cnt>3){
                if (verbose>1){
                    fprintf(stderr, "[debug::%s]     ~ neighbor %.*s hinted more haps but might has another way out\n",
                                __func__, (int)Get_NAME_LENGTH(R_INF, w), Get_NAME(R_INF, w));
                }    
                need_to_protect = 0;
                goto finish;
            }

            if (verbose>1){
                fprintf(stderr, "[debug::%s]     ~ neighbor %.*s hinted more haps\n",
                            __func__, (int)Get_NAME_LENGTH(R_INF, w), Get_NAME(R_INF, w));
            }
            need_to_protect = 1;
            goto finish;
        }else{
            if (verbose>1){
                fprintf(stderr, "[debug::%s]     ~ neighbor %.*s ok\n",
                            __func__, (int)Get_NAME_LENGTH(R_INF, w), Get_NAME(R_INF, w));
            }
        }
        

    }

    
finish:
    if (need_to_protect>0){  // drop arcs between the shorter read and its "parent" reads
        if (verbose){
            fprintf(stderr, "[debug::%s]     protect %.*s ; total %d weak %d\n", __func__, 
                            (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                            nb_total_ovlp, nb_weak_ovlp);
            fprintf(stderr, "[debug::%s]       traceback - the parent reads:\n", __func__);
            for (i=0; i<buf_len; i++){
                fprintf(stderr, "[debug::%s]         %.*s\n", __func__,
                            (int)Get_NAME_LENGTH(R_INF, hits[i]->tn), Get_NAME(R_INF, hits[i]->tn));
            }
        }
        d->counter++;
        for (i=0; i<buf_len; i++){
            hits[i]->del = 1;
            delete_single_edge_both_dir(sources, coverage_cut, (uint32_t)i_r, hits[i]->tn);
            
        }
    }else if (need_to_protect==0){
        if (verbose){
            fprintf(stderr, "[debug::%s]     contained but not protecting %.*s ; total %d weak %d\n", __func__, 
                            (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                            nb_total_ovlp, nb_weak_ovlp);
            fprintf(stderr, "[debug::%s]       traceback - the parent reads:\n", __func__);
            for (i=0; i<buf_len; i++){
                fprintf(stderr, "[debug::%s]         %.*s\n", __func__,
                            (int)Get_NAME_LENGTH(R_INF, hits[i]->tn), Get_NAME(R_INF, hits[i]->tn));
            }
        }
    }else{  // read not contained
        ;
    }

    free(buf);
    free(hits);
}


static void hamt_hit_contained_worker(void *data, long i_r, int tid){  // callback for kt_for()
    // NOTE / might have bugs
    //    Prevent cases from being identified as containment too aggressively will result
    //     in unexpected results. However, need more test set to further improve.
    //    Current solution is to 
    //      1) restrict treatments to reads seemingly from low coeverge regions, 
    //         where it's more likely to see "broken haplotypes" - if we phase it truthfully,
    //         the region ends up fragmented; but if we allow chimera, it can be represented by fewer
    //         contigs or even a single contig. For meta this is preferred, for others, maybe not so much.
    //      2) arbitrarily set thresholds on ovec error counts. The impirical error count distributions
    //         is a light tailed curve, most reads had <50 bases corrected.
    // NOTE about containment
    //    The case
    /*
      read1  ============================
      read2      ============--
    */
   //      ...isn't considered as containment since the right side of read2 doesn't match read1.
    int verbose = 0;  // debug
    dbgmsg_t msg;
    char tmp_msg[500];
    if (verbose){
        hamt_dbgmsg_init(&msg);
    }

    hitcontain_aux_t *d = (hitcontain_aux_t*)data;
    ma_hit_t_alloc *sources = d->sources;
    ma_hit_t_alloc *reverse_sources = d->reverse_sources;
    long long n_read = d->n_read;
    uint64_t *readLen = d->readLen;
    ma_sub_t *coverage_cut = d->coverage_cut;


    // aux
    int nb_reads_contained;
    ma_hit_t *h, *hh;
    kvec_t(float) buf;
    kvec_t(uint32_t) IDs_contained;
    kv_init(buf);
    kv_init(IDs_contained);
    int flag_passed, flag_total;
    float current_error_rate, tmp_error_rate;
    uint16_t current_error_count, tmp_error_count;
    uint32_t tn_flank;
    int readtype;

    int tmp_n = 0, bp_uncovered = 0;
    uint32_t prev_end, prev_start, tmp_start, tmp_end;


    // for interval merge
    kvec_t(uint64_t) spans;
    kv_init(spans);


    buf.n = 0;  // TODO: buffer reset, but this is ugly
    IDs_contained.n = 0;  // TODO: buffer reset, but this is ugly

    if (R_INF.mask_readnorm[i_r] & 1){goto finish;}
    if (verbose){
        sprintf(tmp_msg, "[debug::%s] at read %.*s, homo targets %d\n", 
                            __func__, (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                            (int)sources[i_r].length);
    }

    // topo check
    // readtype = hamt_ovlp_is_end_read(sources, i_r, readLen, 500);  // dont do
    // if (hamt_ovlp_is_contained_read2(sources, i_r, readLen)){continue;}  // dont do
    nb_reads_contained = hamt_ovlp_contains_howmany_reads(sources, i_r, readLen);
    if (nb_reads_contained==0){goto finish;}

    // "coverage" check: do nothing if the region seems well covered (TODO: less arbitrary criteria if more examples found)
    if ((sources[i_r].length/*+reverse_sources[i_r].length*/)>50){
        goto finish;
    }else{
        // (expecting most overalps to be symmetric, not checking from the other direction here)
        ;
    }

    // check number of bases corrected for reads contained by the queryread
    // NOTE: using very crude criteria for now (Dec 2 2020)
    flag_passed = flag_total = 0;
    current_error_count = hamt_count_error(i_r);
    current_error_rate = (float)100*current_error_count/readLen[i_r];
    if (current_error_count<50 || current_error_rate<0.05){goto finish;}
    if (verbose) {
        sprintf(tmp_msg, "[debug::%s] at read %.*s, homo targets %d\n", 
                            __func__, (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r),
                            (int)sources[i_r].length);
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
        sprintf(tmp_msg, "[debug::%s]    current error count is %d, percentage rate %.2f\n", 
                                __func__, (int)current_error_count, current_error_rate);
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
    }
    for (int i=0; i<sources[i_r].length; i++){
        h = &sources[i_r].buffer[i];
        if (h->del) {continue;}
        if ( (h->ts==0 && h->te==readLen[h->tn]) /*|| (h->ts==readLen[h->tn] && h->te==0)*/ ){  // target is contained
            kv_push(uint32_t, IDs_contained, h->tn);
            kv_push(uint64_t, spans, h->qns<<32 | (uint64_t)h->qe );

            tmp_error_count = hamt_count_error(h->tn);
            tmp_error_rate = (float)100*tmp_error_count/readLen[h->tn];
            if (current_error_rate<tmp_error_rate*2){
                if (verbose) {
                    sprintf(tmp_msg, "[debug::%s]    failed an ec check: target %.*s, cnt %d\n", __func__,
                                        (int)Get_NAME_LENGTH(R_INF, h->tn), Get_NAME(R_INF, h->tn), tmp_error_count);
                    hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
                }
                flag_passed++;  // query read isn't too erroneous
            }
            flag_total++;
        }else if (h->ts==0 || h->te==readLen[h->tn]){  // ideal overlap
            kv_push(uint64_t, spans, h->qns<<32 | (uint64_t)h->qe );
        }else{  // non-ideal overlap; TODO: better heuristic if more examples found
            uint32_t adj_s, adj_e;  // "adjusted start/end"
            if (h->rev==0){
                adj_s = ((uint32_t)h->qns - h->ts)<0 ? 0 : ((uint32_t)h->qns - h->ts);
                adj_e = h->qe+(readLen[h->tn]-h->te)>readLen[i_r]? readLen[i_r] : h->qe+(readLen[h->tn]-h->te);
            }else{
                adj_s = (uint32_t)h->qns-(readLen[h->tn]-h->te)<0 ? 0:((uint32_t)h->qns-(readLen[h->tn]-h->te));
                adj_e = h->qe+h->ts>readLen[i_r]? readLen[i_r] : h->qe+h->ts;
            }
            kv_push(uint64_t, spans, ((uint64_t)adj_s)<<32 | (uint64_t)adj_e);
        }
    }


    if ((float)flag_passed/flag_total > 0.5){
        if (verbose) {
            sprintf(tmp_msg, "[debug::%s] query ok won't cut, passed %d, total %d\n", __func__, flag_passed, flag_total);
            hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
        }
        goto finish;
    }
    if (verbose) {
        sprintf(tmp_msg, "[debug::%s] query read erroneous, consider cutting\n", __func__);
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
    }

    // check if the current read could be fully covered by its contained reads and outgoing targets
    // TODO: oneshot, but a wrapper would look nicer
    radix_sort_hamtu64(spans.a, spans.a+spans.n);
    if (spans.n==0){
        sprintf(tmp_msg, "[E::%s] span buffer is empty, this shouldn't happen\n", __func__);
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
        goto finish;
    }
    
    prev_end=(uint32_t)spans.a[0], prev_start=(uint32_t)(spans.a[0]>>32);
    if (spans.n==1){
        tmp_n = 1;
    }else{
        for (int i=1; i<spans.n; i++){
            tmp_start = (uint32_t)(spans.a[i]>>32);
            tmp_end = (uint32_t)spans.a[i];
            if (tmp_start>prev_end){
                bp_uncovered += tmp_start - prev_end;
                // log (recycle the vector's buffer)
                spans.a[tmp_n] = ((uint64_t)prev_start)<<32 | (uint64_t)prev_end;
                prev_end = tmp_end;
                prev_start = tmp_start;
                tmp_n++;
            }else{
                prev_end = tmp_end;
            }
        }
        // last one
        spans.a[tmp_n] = ((uint64_t)prev_start)<<32 | (uint64_t)prev_end;
        tmp_n++;
    }
    if (bp_uncovered>500 || spans.n>=3){
        if (verbose) {
            sprintf(tmp_msg, "[debug::%s]     query erroneous but itself ill covered,don't cut\n", __func__);
            hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
        }
        goto finish;
    }
    
    // if we reach here, it's safe to favor the contained reads over the current longer read
    delete_all_edges(sources, coverage_cut, (uint32_t)i_r);
    if (verbose){
        sprintf(tmp_msg, "[debug::%s] ditched read %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r));
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
        sprintf(tmp_msg, "[debug::%s]     (there were %d contained reads)\n", __func__, (int)IDs_contained.n);
        hamt_dbgmsg_append(&msg, tmp_msg, strlen(tmp_msg));
    }
    d->counter++;

finish:
    kv_destroy(buf);
    kv_destroy(IDs_contained);
    kv_destroy(spans);
    if (verbose){
        if (!hamt_dbgmsg_is_empty(&msg)){
            fprintf(stderr, "%s", msg.a);
        }
        hamt_dbgmsg_destroy(&msg);
    }
}


void hamt_hit_contained_multi(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long n_read, uint64_t *readLen,
                        ma_sub_t *coverage_cut)
{
    // multithread ver of hamt_hit_contained (removed since r025->r026)
    // NOTE:
    //    This alternative function was introduced because sometimes shorter reads might be considered
    //      contained by a longer but actually erroneous read (not sequencing error, but shared subsequences
    //      between similar species etc), which, if ends up in the unitig, might terminate the single walk incorrectly.
    //     First, for each read that contains other reads && is not contained by any other reads,
    //            compare its number of bases corrected to the reads it contains. If the contained reads are
    //            obviously better, *consider* drop the read (see code for additional criteria), or just mark
    //            some overlaps as removed.
    //     Then, perform the regualr contained reads removal.
    double startTime = Get_T();

    hitcontain_aux_t aux;
    aux.sources = sources;
    aux.reverse_sources = reverse_sources;
    aux.n_read = n_read;
    aux.readLen = readLen;
    aux.coverage_cut = coverage_cut;
    aux.counter = 0;

    kt_for(asm_opt.thread_num, hamt_hit_contained_worker, &aux, (long)n_read);
    // kt_for(1, hamt_hit_contained_worker, &aux, (long)n_read);

    fprintf(stderr, "[M::%s] treated roughly %d spots, used %.1fs\n", 
        __func__, aux.counter, Get_T()-startTime);
}

void hamt_hit_contained_drop_singleton_multi(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long n_read, uint64_t *readLen,
                        ma_sub_t *coverage_cut){
    double startTime = Get_T();

    hitcontain_aux_t aux;
    aux.sources = sources;
    aux.reverse_sources = reverse_sources;
    aux.n_read = n_read;
    aux.readLen = readLen;
    aux.coverage_cut = coverage_cut;
    aux.counter = 0;

    // kt_for(asm_opt.thread_num, hamt_hit_contained_drop_singleton_worker, &aux, (long)n_read);
    kt_for(asm_opt.thread_num, hamt_hit_contained_drop_singleton_worker_v2, &aux, (long)n_read);
    // kt_for(1, hamt_hit_contained_drop_singleton_worker_v2, &aux, (long)n_read);

    fprintf(stderr, "[M::%s] treated roughly %d spots, used %.1fs\n", 
                __func__, aux.counter, Get_T()-startTime);


}



#if 0
static void hamt_hit_drop_high_mismatch_arcs_worker(void *data, long i_r, int tid){  // callback for kt_for()
    int verbose = 0;
    hitcontain_aux_t *d = (hitcontain_aux_t*)data;
    ma_hit_t_alloc *sources = d->sources;
    long long n_read = d->n_read;
    uint64_t *readLen = d->readLen;

    ma_hit_t *h;
    float nm;
    UC_Read read1, read2;
    init_UC_Read(&read1);
    init_UC_Read(&read2);

    for (int i=0; i<sources[i_r].length; i++){
        h = &sources[i_r].buffer[i];
        nm = (float)hamt_count_mismatches_in_an_exact_overlap(h, &read1, &read2);
        if (verbose){
            fprintf(stderr, "[debug::%s] qn %.*s (%d) tn %.*s (%d) rev %d nm %d\n", __func__,
                            (int)Get_NAME_LENGTH(R_INF, i_r), Get_NAME(R_INF, i_r), (int)i_r,
                            (int)Get_NAME_LENGTH(R_INF, h->tn), Get_NAME(R_INF, h->tn), (int)h->tn,
                            (int)h->rev,
                             (int)nm);
        }
        if (nm>50 /*|| nm/readLen[i_r]>*/){
            h->del = 1;
            d->counter = d->counter+1;
        }
    }
    destory_UC_Read(&read1);
    destory_UC_Read(&read2);
}
void hamt_hit_drop_high_mismatch_arcs(ma_hit_t_alloc* sources, 
                                      long long n_read, uint64_t *readLen){
    hitcontain_aux_t aux;
    aux.sources = sources;
    aux.n_read = n_read;
    aux.readLen = readLen;
    aux.counter = 0;

    kt_for(asm_opt.thread_num, hamt_hit_drop_high_mismatch_arcs_worker, &aux, (long)n_read);
    fprintf(stderr, "[M::%s] dropped %d arcs\n", __func__, aux.counter);
}
#endif



void ma_hit_cut(ma_hit_t_alloc* sources, long long n_read, uint64_t* readLen, 
long long mini_overlap_length, ma_sub_t** coverage_cut)
{
    double startTime = Get_T();
    int nb_treatedA = 0, nb_treatedB=0;
	size_t i, j;
    ma_hit_t* p;
    ma_sub_t* rq;
    ma_sub_t* rt;
    long long rLen = 0;
    for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        rLen = 0;
        for (j = 0; j < sources[i].length; j++)
        {
            ///this is a overlap
            p = &(sources[i].buffer[j]);
            if(p->del) continue;

            rq = &((*coverage_cut)[Get_qn(*p)]);
            rt = &((*coverage_cut)[Get_tn(*p)]);
            ///if any of target read and the query read has no enough coverage
		    if (rq->del || rt->del) continue;
            int qs, qe, ts, te;




            ///target and query in different strand
            if (p->rev) 
            {
                /**
                here is an example in different strand:
                                      
                                       (te) (rt->e)           (rt->s)    (ts)
                                        |      |                 |        |
                target ----------------------------------------------------
                                        ------------------------------------------- query
                                        qs                                qe
                **/
                qs = p->te < rt->e? Get_qs(*p): Get_qs(*p) + (p->te - rt->e);
                qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
                ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
                te = Get_qs(*p) > rq->s? p->te : p->te - (rq->s - Get_qs(*p));
            }
            else ///target and query in same strand
            { 
                /**
                note: ts is the targe start in this overlap, 
                while rt->s is the high coverage start in the whole target (not only in this overlap)
                so this line is to normalize the qs in quey to high coverage region
                **/
                //(rt->s - p->ts) is the offset        
                qs = p->ts > rt->s? Get_qs(*p): Get_qs(*p) + (rt->s - p->ts); 
                //(p->te - rt->e) is the offset
                qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
                //(rq->s - Get_qs(*p) is the offset
                ts = Get_qs(*p) > rq->s? p->ts : p->ts + (rq->s - Get_qs(*p));
                //(p->qe - rq->e) is the offset
                te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
		    }

            

            //cut by self coverage
            //and normalize the qs, qe, ts, te by rq->s and rt->e
            qs = ((uint32_t)qs > rq->s? qs : rq->s) - rq->s;
            qe = ((uint32_t)qe < rq->e? qe : rq->e) - rq->s;
            ts = ((uint32_t)ts > rt->s? ts : rt->s) - rt->s;
            te = ((uint32_t)te < rt->e? te : rt->e) - rt->s;

            if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length) 
            {
                ///p->qns = p->qns>>32<<32 | qs;
                p->qns = p->qns>>32;
                p->qns = p->qns << 32;
                p->qns = p->qns | qs;

                p->qe = qe;
                p->ts = ts;
                p->te = te;
                p->del = 0;
                rLen++;
		    }
            else
            {
                p->del = 1;
                delete_single_edge(sources, (*coverage_cut), Get_tn(*p), Get_qn(*p));
                nb_treatedA++;
                ///delete_all_edges(paf, coverage_cut, i); 
            }
        }

        if(rLen == 0)
        {
            (*coverage_cut)[i].del = 1;
            nb_treatedB++;
        }
    }

    fprintf(stderr, "[M::%s] typeA %d , typeB %d, used %.1f\n", __func__, nb_treatedA, nb_treatedB,  Get_T()-startTime);
    
}


// typedef struct{
//     ma_hit_t_alloc *sources;
//     ma_hit_t_alloc *reverse_sources;
//     long long n_read;
//     uint64_t *readLen;
//     int cov_nbreads;
//     int counter;  // might RACE but doesn't matter
// }hitres_hapgap_t;  // "hit rescue haplotype gap"




/**********************************
 * Filter short potential unitigs *
 **********************************/
static inline int asg_is_utg_end(const asg_t *g, uint32_t v, uint64_t *lw)
{
	uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
	int i, i0 = -1;
	asg_arc_t *aw, *av = asg_arc_a(g, v^1);

	for (i = nv = 0; i < (int)nv0; ++i)  // count still existing in-arcs
		if (!av[i].del) i0 = i, ++nv;

	///end without any out-degree
	if (nv == 0) return ASG_ET_TIP; // tip
	if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs
	if (lw) *lw = av[i0].ul<<32 | av[i0].v;
	w = av[i0].v ^ 1;
	nw0 = asg_arc_n(g, w);
	aw = asg_arc_a(g, w);
	for (i = nw = 0; i < (int)nw0; ++i)
		if (!aw[i].del) ++nw;
	if (nw != 1) return ASG_ET_MULTI_NEI;
	return ASG_ET_MERGEABLE;
}


int asg_extend(const asg_t *g, uint32_t v, int max_ext, asg64_v *a)
{
	int ret;
	uint64_t lw;
	a->n = 0;
	kv_push(uint64_t, *a, v);
	do {
        /**
		 test (v) and (v^1), 
         if the out-degrees of both (v) and (v^1) are 1, ret == 0
        **/
		ret = asg_is_utg_end(g, v^1, &lw);
        /**
        #define ASG_ET_MERGEABLE 0
        #define ASG_ET_TIP       1
        #define ASG_ET_MULTI_OUT 2
        #define ASG_ET_MULTI_NEI 3
        **/
		if (ret != 0) break;
		kv_push(uint64_t, *a, lw);
		/**
		 ret == 0 means: 
		   v^1 and is the only prefix of (uint32_t)lw, 
		   and (uint32_t)lw is the only prefix of v^1
		**/
		
		v = (uint32_t)lw;
	} while (--max_ext > 0);
	return ret;
}


static inline int asg_is_single_edge(const asg_t *g, uint32_t v, uint32_t start_node)
{// hamt note: function name does not describe what it'd do?
 //        (count parent nodes (regardless of deletion state))
 //        (then why not use asg_arc_n(g, v^1) directly???)

	uint32_t nv, nv0 = asg_arc_n(g, v^1);
	int i;
	asg_arc_t *av = asg_arc_a(g, v^1);

    int flag = 0;
	///if this arc has not been deleted
	for (i = nv = 0; i < (int)nv0; ++i)
    {
        ///if (!av[i].del)
        {
            ++nv;
            if(av[i].v>>1 == start_node)  // (hamt note)i think this is always true unless graph is broken?
            {
                flag = 1;
            }
        } 
    }
		
    if(flag == 0)
    {
        fprintf(stderr, "****ERROR\n");
    }

    return nv;
}

void debug_info_of_specfic_node(char* name, asg_t *g, R_to_U* ruIndex, char* command)
{
    fprintf(stderr, "\n\n\n");
    uint32_t v, n_vtx = g->n_seq * 2, queryLen = strlen(name), flag = 0, contain_rId, is_Unitig;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(queryLen == Get_NAME_LENGTH(R_INF, (v>>1)) && memcmp(name, Get_NAME(R_INF, (v>>1)), Get_NAME_LENGTH(R_INF, (v>>1))) == 0)
        {
            if(flag == 0) fprintf(stderr, "\nafter %s\n", command);
            fprintf(stderr, "****************graph ref_read: %.*s, id: %u, dir: %u****************\n", 
            (int)Get_NAME_LENGTH(R_INF, (v>>1)), Get_NAME(R_INF, (v>>1)), v>>1, v&1);
            if(g->seq[v>>1].del)
            {
                get_R_to_U(ruIndex, (v>>1), &contain_rId, &is_Unitig);
                if(contain_rId != (uint32_t)-1 && is_Unitig != 1)
                {
                    fprintf(stderr, "read is deleted as a contained read by: %.*s\n contain_rId: %u, del: %u\n", 
                    (int)Get_NAME_LENGTH(R_INF, contain_rId), Get_NAME(R_INF, contain_rId), contain_rId, g->seq[contain_rId].del);
                }
                else
                {
                    fprintf(stderr, "read has already been deleted.\n");
                }
                
                return;
            }

            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t i, nv = asg_arc_n(g, v);
            for (i = 0; i < nv; ++i)
            {
                fprintf(stderr, "target: %.*s, el: %u, strong: %u, ol: %u, del: %u\n", 
                (int)Get_NAME_LENGTH(R_INF, (av[i].v>>1)), 
                Get_NAME(R_INF, (av[i].v>>1)),
                av[i].el, av[i].strong, av[i].ol, av[i].del);
            }
            flag = 1;
        }
    }
}


asg_t *ma_sg_gen(const ma_hit_t_alloc* sources, long long n_read, const ma_sub_t *coverage_cut, 
int max_hang, int min_ovlp)
{
    double startTime = Get_T();
	size_t i, j;
	asg_t *g;
	///just calloc
	g = asg_init();

	///add seq to graph, seq just save the length of each read
	for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        ///if a read has been deleted, should we still add them?
		asg_seq_set(g, i, coverage_cut[i].e - coverage_cut[i].s, coverage_cut[i].del);
        g->seq[i].c = coverage_cut[i].c;  // i think it should be PRIMARY_LABEL i.e. 0?
	}

    g->seq_vis = (uint8_t*)calloc(g->n_seq*2, sizeof(uint8_t));

    for (i = 0; i < (uint64_t)n_read; ++i) 
    {
        for (j = 0; j < sources[i].length; j++)
        {
            int r;
		    asg_arc_t t, *p;
		    const ma_hit_t *h = &(sources[i].buffer[j]);
            if(h->del) continue;

            //high coverage region [sub[qn].e, sub[qn].s) in query
            int ql = coverage_cut[Get_qn(*h)].e - coverage_cut[Get_qn(*h)].s;
            //high coverage region [sub[qn].e, sub[qn].s) in target
            int tl = coverage_cut[Get_tn(*h)].e - coverage_cut[Get_tn(*h)].s;
		    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
            /**
            #define MA_HT_INT       (-1)
            #define MA_HT_QCONT      (-2)
            #define MA_HT_TCONT      (-3)
            #define MA_HT_SHORT_OVLP (-4)
            the short overlaps and the overlaps with contain reads have already been removed
            here we should have overhang
            so r should always >= 0
            **/
            if (r >= 0) 
            {
                ///push node?
                p = asg_arc_pushp(g);
                *p = t;
		    }
            else
            {
                fprintf(stderr, "[E::%s]error, dangling bad overlap. r is %d , ql %d tl %d\n", 
                        __func__, r, ql, tl);
                ma_hit_t *h_san=0;
                uint32_t tn = Get_tn(*h);
                uint32_t qn = Get_qn(*h);
                for (int jj=0; jj<sources[tn].length; jj++){
                    if (Get_tn(sources[tn].buffer[jj])==qn) {
                        h_san = &sources[tn].buffer[jj];
                        break;
                    }
                }
                if (!h_san){
                    fprintf(stderr, "sancheck: self is %d, target %d, the other way around is not found!\n", 
                            (int)qn, (int)tn);
                }else{
                    fprintf(stderr, "(sancheck: self is %d, target %d, r of the other way around is %d , self coverage cut del=%d, opposite covcut del=%d)\n", 
                        (int)qn, (int)tn, 
                        ma_hit2arc(h_san, tl, ql, max_hang, asm_opt.max_hang_rate, min_ovlp, &t),
                        (int)coverage_cut[Get_qn(*h)].del, (int)coverage_cut[Get_tn(*h)].del
                       );
                }
            }
        }
    }

    asg_cleanup(g);
    g->r_seq = g->n_seq;

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] used %.1f s\n", __func__, Get_T()-startTime);
    }

	return g;
}




///pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
///note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_finder_with_del_advance(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	if (g->seq[v0>>1].del) return 0;
	if ((uint32_t)g->idx[v0] < 2) return 0; // effectively asg_arc_n(n0)
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S);
        uint32_t d = b->a[v].d;  // "the shortest distance from the initial vertex"
        uint32_t c = b->a[v].c;  // "max count of positive reads"
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
        /****************************may have bugs********************************/
		///assert(nv > 0);
        /****************************may have bugs********************************/

		///all out-edges of v
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul;
			binfo_t *t = &b->a[w];  // vertex info of w
			///that means there is a circle, directly terminate the whole bubble poping
			if (w == v0)  // is looping, abort
            {
                //fprintf(stderr, "n_pop error1\n");
                goto pop_reset;
            } 
            
			///if this edge has been deleted
            /****************************may have bugs********************************/
			///if (av[i].del) continue;
            /****************************may have bugs********************************/

			///push the edge // b->e is "visited edges/arcs"
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist)
            {
                //fprintf(stderr, "n_pop error2\n");
                break; // too far
            } 
            
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l, t->c = c + 1;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
                // my note: swap parent node & update numbers if the current path is longer than the record
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            ///fprintf(stderr, "n_pop error3\n");
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
    ///fprintf(stderr, "n_pop: %d\n", n_pop);
	return n_pop;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_finder_without_del_advance(asg_t *g, uint32_t v0, int max_dist, 
buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
    if(count_out_without_del(g, v0) < 2) return 0; // no bubbles


	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
        /****************************may have bugs********************************/
		///assert(nv > 0);
        /****************************may have bugs********************************/

		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			if (w == v0)
            {
                //fprintf(stderr, "n_pop error1\n");
                goto pop_reset;
            } 
            
			///if this edge has been deleted
            /****************************may have bugs********************************/
			if (av[i].del) continue;
            /****************************may have bugs********************************/

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist)
            {
                //fprintf(stderr, "n_pop error2\n");
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l, t->c = c + 1;
				///incoming edges of w
				t->r = count_out_without_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            ///fprintf(stderr, "n_pop error3\n");
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
    ///fprintf(stderr, "n_pop: %d\n", n_pop);
	return n_pop;
}


// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
//note!!!!!!!! here we don't exculde the deleted edges
static uint64_t asg_bub_end_finder_with_del_advance(asg_t *g, uint32_t* v_Ns, uint32_t occ, 
int max_dist, buf_t *b, uint32_t exculde_init, uint32_t exclude_node, uint32_t* sink)
{
	uint32_t i, j, n_pending = 0;
	uint64_t n_pop = 0;
    ///S saves nodes with all incoming edges visited
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    for (j = 0; j < occ; j++)
    {
        ///if this node has been deleted
        if (g->seq[v_Ns[j]>>1].del) return 0; // already deleted
        ///for each node, b->a saves all related information
        b->a[v_Ns[j]].c = b->a[v_Ns[j]].d = 0;
        ///b->S is the nodes with all incoming edges visited
        kv_push(uint32_t, b->S, (v_Ns[j]<<1)|exculde_init);
    }
    

	
	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), f = v & (uint32_t)1;
        v = v >> 1;
        uint32_t d = b->a[v].d, c = b->a[v].c;


		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		
		///all out-edges of v
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul;
			binfo_t *t = &b->a[w];

            for (j = 0; j < occ; j++)
            {
                if(w == v_Ns[j]) goto pop_reset;
            }
            


            if(f && (exclude_node) == (w)) continue;

            ///if (av[i].del) continue;

			///push the edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);

			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist)
            {
                break; // too far
            } 
            


			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p means the in-node of w is v
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l, t->c = c + 1;
				///incoming edges of w
				t->r = count_out_with_del(g, w^1);
                ///t->r = count_out_without_del(g, w^1);
				++n_pending;
			} else { // visited before
				///c seems the max weight of node
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
				if (d + l < t->d) t->d = d + l; // update dist
			}
            /****************************may have bugs********************************/
			///assert(t->r > 0);
            /****************************may have bugs********************************/
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				//if (x) kv_push(uint32_t, b->S, w);
                if (x) kv_push(uint32_t, b->S, w<<1);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0)
        {
            goto pop_reset;
        } 
        
	} while (b->S.n > 1 || n_pending);

    (*sink) = b->S.a[0]>>1;
	///asg_bub_backtrack(g, v0, b);
	///n_pop = 1 | (uint64_t)b->T.n<<32;
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}




int if_node_exist(uint32_t* nodes, uint32_t length, uint32_t query)
{
    uint32_t i;
    for (i = 0; i < length; ++i)
    {
        if((nodes[i]>>1) == query)
        {
            return 1;
        }
    }

    return 0;
}




long long single_edge_length(asg_t *g, uint32_t begNode, uint32_t endNode, long long edgeLen)
{
    
    uint32_t v = begNode;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    long long rLen = 0;

    while (rLen < edgeLen && nv == 1)
    {
        rLen++;
        if((av[0].v>>1) == endNode)
        {
            return rLen;
        }

        if(asg_is_single_edge(g, av[0].v, v>>1) != 1)
        {
            return -1;
        }

        v = av[0].v;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);      
    }

    return -1;   

}

uint32_t detect_single_path(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b)
{
    
    uint32_t v = begNode;
    uint32_t nv, rnv;
    asg_arc_t *av;
    (*Len) = 0;


    while (1)
    {
        (*Len)++;
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);   
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(nv == 0)
        {
            return END_TIPS;
        }

        if(nv == 2)
        {
            return TWO_OUTPUT;
        }

        if(nv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, nv=1
        ///rnv must >= 1
        rnv = asg_is_single_edge(g, av[0].v, v>>1);
        v = av[0].v;
        (*endNode) = v;
        if(rnv == 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(rnv > 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

int detect_bubble_end(asg_t *g, uint32_t begNode1, uint32_t begNode2, uint32_t* endNode, 
long long* minLen, buf_t* b)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(detect_single_path(g, begNode1, &e1, &l1, b) == TWO_INPUT 
       && 
       detect_single_path(g, begNode2, &e2, &l2, b) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            return 1;
        }
    }

    return 0;
}


int detect_simple_bubble(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* minLen, buf_t* b)
{  // the simplest and very strict bubble 
   //(does not allow non-single edge along the buble edges, does not allow multi out/in at the begin/end node)
    uint32_t e1, e2;
    long long l1, l2;

    if(asg_arc_n(g, begNode) != 2)
    {
        return 0;
    }

    if(asg_is_single_edge(g, asg_arc_a(g, begNode)[0].v, begNode>>1)!=1
      || 
       asg_is_single_edge(g, asg_arc_a(g, begNode)[1].v, begNode>>1)!=1)
    {
        return 0;
    }

    if(b) kv_push(uint32_t, b->b, begNode>>1);

    if(detect_single_path(g, asg_arc_a(g, begNode)[0].v, &e1, &l1, b) == TWO_INPUT 
       && 
       detect_single_path(g, asg_arc_a(g, begNode)[1].v, &e2, &l2, b) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            (*minLen)++;
            return 1;
        }
    }

    return 0;
}


uint32_t detect_single_path_with_single_bubbles(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* Len, buf_t* b, uint32_t max_ext)
{   // detect single path, but allow simple&strict bubble(s)
    uint32_t v = begNode;
    uint32_t nv, rnv;
    asg_arc_t *av;
    long long bLen;
    long long pre_b_n = 0;
    (*Len) = 0;
  
    while (1)
    {
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);   
        (*endNode) = v;
        (*Len)++;

        if((*Len) > max_ext)
        {
            return LONG_TIPS_UNDER_MAX_EXT;
        }

        if(b) kv_push(uint32_t, b->b, v>>1);  // push to visited vertices

        if(nv == 0)
        {
            return END_TIPS;
        }

        if(nv == 2)
        {
            if(b) pre_b_n = b->b.n;
            // the bifurcation isn't the start of a simple&strict bubble
            if(!detect_simple_bubble(g, v, &v, &bLen, b)) 
            {
                if(b) b->b.n = pre_b_n;
                return TWO_OUTPUT;
            }
            // bifurcation is the start of a simple&strict bubble
            // it's ok to continue
            // v is updated to the end of the bubble
            (*Len) = (*Len) + bLen - 2;
            continue;
        }

        if(nv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, nv=1
        ///rnv must >= 1
        rnv = asg_is_single_edge(g, av[0].v, v>>1);
        v = av[0].v;
        (*endNode) = v;
        if(rnv == 2)
        {
            if(b) kv_push(uint32_t, b->b, v>>1);
            (*Len)++;
            return TWO_INPUT;
        }

        if(rnv > 2)
        {
            if(b) kv_push(uint32_t, b->b, v>>1);
            (*Len)++;
            return MUL_INPUT;
        } 

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        }  
    }

    return LONG_TIPS;   
}

int detect_bubble_end_with_bubbles(asg_t *g, uint32_t begNode1, uint32_t begNode2, 
uint32_t* endNode, long long* minLen, buf_t* b)
{  // begNode1 and begNode2 meet at the same vertex after a single path (which allows simple&strict bubbles) of their own downstream traversal
    uint32_t e1, e2;
    long long l1, l2;
    if(detect_single_path_with_single_bubbles(g, begNode1, &e1, &l1, b, (uint32_t)-1) == TWO_INPUT 
       && 
       detect_single_path_with_single_bubbles(g, begNode2, &e2, &l2, b, (uint32_t)-1) == TWO_INPUT)
    {
        if(e1 == e2)
        {
            (*endNode) = e1;
            (*minLen) = (l1 <= l2)? l1: l2;
            return 1;
        }
    }
    return 0;
}


int detect_mul_bubble_end_with_bubbles(asg_t *g, uint32_t* begs, uint32_t occ, 
uint32_t* endNode, long long* minLen, buf_t* b)
{
    uint32_t e, flag, e_s;
    long long l, i, l_s;

    if(occ < 1) return 0;

    flag = detect_single_path_with_single_bubbles(g, begs[0], &e, &l, b, (uint32_t)-1);

    if(flag == TWO_INPUT || flag == MUL_INPUT)
    {
        e_s = e;
        l_s = l;
    }
    else
    {
        return 0;
    }
    


    for (i = 1; i < occ; i++)
    {
        flag = detect_single_path_with_single_bubbles(g, begs[i], &e, &l, b, (uint32_t)-1);
        if(flag == TWO_INPUT || flag == MUL_INPUT)
        {
            if(e != e_s) return 0;
            if(l < l_s) l_s = l;
        }
        else
        {
            return 0;
        }
    }

    if(asg_arc_n(g, e_s^1) == occ)
    {
        (*endNode) = e_s;
        (*minLen) = l_s;
        return 1;
    }
    
    return 0;
}


int detect_bubble_with_bubbles(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* minLen, 
buf_t* b, uint32_t max_ext)
{
    uint32_t e1, e2;
    long long l1, l2;

    if(asg_arc_n(g, begNode) != 2)
    {
        return 0;
    }

    if(asg_is_single_edge(g, asg_arc_a(g, begNode)[0].v, begNode>>1)!=1
      || 
       asg_is_single_edge(g, asg_arc_a(g, begNode)[1].v, begNode>>1)!=1)
    {
        return 0;
    }

    if(b) kv_push(uint32_t, b->b, begNode>>1);


    if(detect_single_path_with_single_bubbles(g, asg_arc_a(g, begNode)[0].v, &e1, &l1, b, max_ext) == TWO_INPUT)
    {
        b->b.n--;
        if(detect_single_path_with_single_bubbles(g, asg_arc_a(g, begNode)[1].v, &e2, &l2, b, max_ext) == TWO_INPUT)
        {
            b->b.n--;
            if(e1 == e2)
            {
                (*endNode) = e1;
                (*minLen) = (l1 <= l2)? l1: l2;
                (*minLen)++;
                return 1;
            }
        }
    }

    return 0;
}

/* NOTE
    Written during r43->r44, dropped because r44 was worse, and that this part doesn't matter...
    The case was: some overlaps were actually exact, but perhaps because seeding problems,
       they were reported to be inexact, only off by 1 or a few base pairs.
       For example, readA and readB overlap by spans [0, a] and [0, b],
       the mistaken case is to report the spans as, say, [0, a+1] and [0, b],
       which will result the pair to be flagged as inexact by ha's if_exact_match.
*/
int hamt_check_exact_match_for_some_cases_worker(ma_hit_t *h, uint64_t *readLen,
                                          UC_Read* q_read, UC_Read* t_read,
                                          int threshold){
	// NOTE
	//     This function should be used in addition the very initial sanchecks/cleanings.
	// RETURN
	//     0 if defintely ok
	//     INT if could match but *the ends* are off by INT bp
	//     -1 if defintely no, given the threshold
	//     -2 pass, no comment (if the arc is deleted or something)
	if (h->del) return -2;
    // fprintf(stderr, "check %d\n", 0);
    if ((h->qns>>32)!=0 && h->ts!=0) return -2;  // no comment for non-dovetail/containment ones

	int l_diff = 0;
	uint32_t ol1, ol2;
	uint32_t qn, tn, qe, te, qs, ts;
	qn = (uint32_t) (h->qns>>32);
	tn = h->tn;
	ol1 = h->qe - (uint32_t)h->qns;
	ol2 = h->te - h->ts;
    qs = (uint32_t)h->qns;
	qe = h->qe;
    ts = h->ts;
	te = h->te;
	char *q_seq;
    char *t_seq;

    if (ol1==ol2) return -2;
    // fprintf(stderr, "check %d\n", 1);
    if (qs!=0 && qe!=readLen[qn]) return -2;
    // fprintf(stderr, "check %d\n", 2);
    if (ts!=0 && te!=readLen[tn]) return -2;
    // fprintf(stderr, "check %d\n", 3);

	if (ol1>ol2){
		l_diff = ol1-ol2;
		if (qs==0) qe-=l_diff;
        else qs+=l_diff;
        ol1 = ol2;
	}else if (ol2>ol1){
		l_diff = ol2-ol1;
		if (ts==0) te-=l_diff;
        else ts+=l_diff;
        ol2 = ol1;
	}
    if (l_diff>threshold) return -1;
    // fprintf(stderr, "check %d\n", 4);

    recover_UC_sub_Read(q_read, qs, ol1, 0, &R_INF, qn);
    recover_UC_sub_Read(t_read, ts, ol1, h->rev, &R_INF, tn);

    if (memcmp(q_read->seq, t_read->seq, ol1)==0){
        // fprintf(stderr, "check %d\n", 5);
        return l_diff;
    }else{
        return -1;
    }
}
int hamt_rescue_inexact_by_1bp(ma_hit_t_alloc* sources, uint64_t *readLen){
    // NOTE
    //     This function only alters the .el marking, won't change the actual alignment entry.
    int verbose = 0;
    UC_Read read1, read2;
    init_UC_Read(&read1);
    init_UC_Read(&read2);
    ma_hit_t *h;
    uint32_t w;
    int status;
    int ret = 0;
    int abnormal = 0;
    for (uint32_t v=0; v<R_INF.total_reads; v++){
        if (verbose) fprintf(stderr, "[debug::%s] at %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, v), Get_NAME(R_INF, v));
        for (int i=0; i<sources[v].length; i++){
            h = &sources[v].buffer[i];
            w = h->tn;
            if (sources[v].buffer[i].el) continue;   // skip overlaps already marked as exact
            if (verbose) fprintf(stderr, "[debug::%s]   target %.*s\n", __func__, 
                                            (int)Get_NAME_LENGTH(R_INF, sources[v].buffer[i].tn), Get_NAME(R_INF, sources[v].buffer[i].tn));
            status = hamt_check_exact_match_for_some_cases_worker(&sources[v].buffer[i], readLen, &read1, &read2, 1);
            if (verbose) fprintf(stderr, "[debug::%s]     > ret %d\n", __func__, status);
            if (status>=0){
                ret++;
                sources[v].buffer[i].el = 1;
                h = get_specific_overlap_handle(sources, w, v);
                if (!h){
                    abnormal++;
                    continue;
                }
                h->el = 1;
            }
        }
    }
    fprintf(stderr, "[M::%s] converted %d arcs (abnormal encountered: %d)\n", __func__, ret, abnormal);
    return ret;
}




int test_triangular_exact(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode, int max_dist, buf_t* bub)
{
    
    uint32_t i, v, w;
    ///int flag0, flag1, node;
    int n_reduced = 0, todel;
    long long NodeLen_first[3];
    long long NodeLen_second[3];

    uint32_t Ns_first[3];
    uint32_t Ns_second[3];
    
    
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        if(av[0].v == av[1].v)
        {
            continue;
        }

        // probe parent status for the two vertices
        NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
        if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
        {
            NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
            NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
        }
        // (require that one of them has 1 parent and the other has 2)
        if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
        {
            continue;
        }
       
        // check the node who has 2 parents
        if(av[NodeLen_first[2]].del == 1) {continue;}
        w = av[NodeLen_first[2]].v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        if(nw != 2)
        {
            fprintf(stderr, "[E::%s]error\n", __func__);
        }
        NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
        if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
        {
            NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
            NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
        }
        // require that one node has one out-edge, another node has two out-edges
        if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
        {
            continue;
        }


        if(if_node_exist(nodes, length, (w>>1)) && ((w>>1) != (endNode>>1)))  // if w>>1 exists in the given list of vertices, and w>>1 is not the endNode
        {
               uint32_t convex1 = 0, convex2 = 0, f1, f2;
               long long l1 = 0, l2 = 0;
               todel = 0;
               f1 = detect_bubble_end_with_bubbles(g, av[0].v, av[1].v, &convex1, &l1, NULL);  // 1 if yes, else 0 
               f2 = detect_bubble_end_with_bubbles(g, aw[0].v, aw[1].v, &convex2, &l2, NULL);
               if(f1 && f2)
               {
                    if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                     ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                     {
                         if(l1 <= min_thres || l2 <= min_thres)
                         {
                             continue;
                         } 

                         todel = 1;
                     } 
               }
               else if(f1)
               {
                  if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                   {
                       if(l1 <= min_thres)
                       {
                           continue;
                       }

                       todel = 1;
                   }   
               }
               else if(f2)
               {
                   if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                   {
                       if(l2 <= min_thres)
                       {
                           continue;
                       }

                       todel = 1;
                   }   
               }

            
               if(todel == 0)
               {
                   if(!f1)
                   {
                       Ns_first[0] = av[0].v; Ns_first[1] = av[1].v;
                       f1 = asg_bub_end_finder_with_del_advance(g, Ns_first, 2, max_dist, 
                       bub, 0, (u_int32_t)-1, &convex1);
                       l1 = min_thres + 10;
                   }

                   if(!f2)
                   {
                       Ns_second[0] = aw[0].v; Ns_second[1] = aw[1].v;
                       f2 = asg_bub_end_finder_with_del_advance(g, Ns_second, 2, max_dist, 
                       bub, 0, (u_int32_t)-1, &convex2);
                       l2 = min_thres + 10;
                   }




                    
                    if(f1 && f2)
                    {
                        if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                        ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                        {
                            if(l1 <= min_thres || l2 <= min_thres)
                            {
                                continue;
                            } 

                            todel = 1;
                        } 
                    }
                    else if(f1)
                    {
                        if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                        {
                            if(l1 <= min_thres)
                            {
                                continue;
                            }

                            todel = 1;
                        }   
                    }
                    else if(f2)
                    {
                        if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                        {
                            if(l2 <= min_thres)
                            {
                                continue;
                            }

                            todel = 1;
                        }   
                    }
                    
               }
               
               
               if(todel)
               {
                   av[NodeLen_first[2]].del = 1;
                    ///remove the reverse direction
                    asg_arc_del(g, av[NodeLen_first[2]].v^1, av[NodeLen_first[2]].ul>>32^1, 1);
                    n_reduced++;
               }
 
        }

    }

    return n_reduced;
}





int find_single_link(asg_t *g, uint32_t link_beg, int linkLen, uint32_t* link_end)
{
    uint32_t v, w;
    v = link_beg^1;
    uint32_t nv, nw;
    asg_arc_t *av;
    int edgeLen = 0;

    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    if(nv != 1)
    {
        return 0;
    }
    v = av[0].v;

    while (edgeLen < linkLen)
    {
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        if(nv != 1)
        {
            return 0;
        }

        w = v^1;
        nw = asg_arc_n(g, w);
        if(nw == 2)
        {
            (*link_end) = w;
            return 1;
        }
        else if(nw > 2)
        {
            return 0;
        }
        
        v = av[0].v;
        edgeLen++;
    }


    return 0;    
}


int if_edge_exist(asg_arc_t* edges, uint32_t length, uint32_t query)
{
    uint32_t i;
    for (i = 0; i < length; ++i)
    {
        if((edges[i].v>>1) == query)
        {
            return 1;
        }
    }

    return 0;
}

int test_quadangular_with_addition_node(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t addition_node_length)
{
    
    uint32_t i, v, w;
    int flag, occ_v_0, occ_v_1, occ_w_0, occ_w_1;
    int n_reduced = 0;
    uint32_t v_out2_node, w_out2_node;
    uint32_t cut_edge_v, cut_edge_w;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 1)
        {
            continue;
        }
        /**********************test first node************************/
        flag = asg_is_single_edge(g, av[0].v, v>>1);
        if(flag != 2)
        {
            continue;
        }
       /**********************test first node************************/

        if(!find_single_link(g, v, addition_node_length, &w))
        {
            continue;
        }
        v = av[0].v^1;
        ///up to now, v and w is the node what we want
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        nv = asg_arc_n(g, v);
        av = asg_arc_a(g, v);
        if(nv!=2 || nw != 2)
        {
            fprintf(stderr, "[E::%s]error\n", __func__);
        }

        if(!if_node_exist(nodes, length, (v>>1)))
        {
            continue;
        }
        if(!if_node_exist(nodes, length, (w>>1)))
        {
            continue;
        }
        /**********************for v************************/
        occ_v_0 = asg_is_single_edge(g, av[0].v, v>>1);
        occ_v_1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(occ_v_0 == occ_v_1)
        {
            continue;
        }
        if(occ_v_0 < 1 || occ_v_0 > 2)
        {
            continue;
        }
        if(occ_v_1 < 1 || occ_v_1 > 2)
        {
            continue;
        }

        if(occ_v_0 == 2)
        {
            v_out2_node = av[0].v^1;
            cut_edge_v = 0;
        }
        else
        {
            v_out2_node = av[1].v^1;
            cut_edge_v = 1;
        }
        if(!if_node_exist(nodes, length, (v_out2_node>>1)))
        {
            continue;
        }
        /**********************for v************************/

        /**********************for w************************/
        occ_w_0 = asg_is_single_edge(g, aw[0].v, w>>1);
        occ_w_1 = asg_is_single_edge(g, aw[1].v, w>>1);
        if(occ_w_0 == occ_w_1)
        {
            continue;
        }
        if(occ_w_0 < 1 || occ_w_0 > 2)
        {
            continue;
        }
        if(occ_w_1 < 1 || occ_w_1 > 2)
        {
            continue;
        }

        if(occ_w_0 == 2)
        {
            w_out2_node = aw[0].v^1;
            cut_edge_w = 0;
        }
        else
        {
            w_out2_node = aw[1].v^1;
            cut_edge_w = 1;
        }
        if(!if_node_exist(nodes, length, (w_out2_node>>1)))
        {
            continue;
        }
        /**********************for w************************/

        
        if(!if_edge_exist(asg_arc_a(g, w_out2_node), asg_arc_n(g, w_out2_node), (v_out2_node>>1)))
        {
            continue;
        }

        if(!if_edge_exist(asg_arc_a(g, v_out2_node), asg_arc_n(g, v_out2_node), (w_out2_node>>1)))
        {
            continue;
        }


        av[cut_edge_v].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, av[cut_edge_v].v^1, av[cut_edge_v].ul>>32^1, 1);

        aw[cut_edge_w].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, aw[cut_edge_w].v^1, aw[cut_edge_w].ul>>32^1, 1);


        n_reduced++;
    }

    return n_reduced;
}

int test_triangular_addition_exact(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode, int max_dist, buf_t* bub)
{
    
    uint32_t i, j, v, w;
    ///int flag0, flag1, node;
    int n_reduced = 0, todel;
    uint32_t Nodes1[2]={0};
    uint32_t Nodes2[2]={0};

    uint32_t Ns_first[2]={0};
    uint32_t Ns_second[2]={0};


    
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 1)
        {
            continue;
        }
        if(asg_is_single_edge(g, av[0].v, v>>1) != 2)
        {
            continue;
        }


        w = v^1;
        asg_arc_t *aw = asg_arc_a(g, w);
        uint32_t nw = asg_arc_n(g, w);
        if(nw != 1)
        {
            continue;
        }
        if(asg_is_single_edge(g, aw[0].v, w>>1) != 2)
        {
            continue;
        }

        if((av[0].v>>1) == (aw[0].v>>1))
        {
            continue;
        }



        Nodes1[0] = av[0].v^1;
        Nodes2[0] = aw[0].v^1;

        
        for(j = 0; j < 2; j++)
        {
            if((asg_arc_a(g, Nodes1[0])[j].v>>1)!= (v>>1))
            {
                Nodes1[1] = asg_arc_a(g, Nodes1[0])[j].v^1;
            }
        }

        for(j = 0; j < 2; j++)
        {
            if((asg_arc_a(g, Nodes2[0])[j].v>>1)!= (v>>1))
            {
                Nodes2[1] = asg_arc_a(g, Nodes2[0])[j].v^1;
            }
        }
        
        if(asg_arc_n(g, Nodes1[1]) != 1 || asg_arc_n(g, Nodes2[1]) != 1)
        {
            continue;
        }

        if((Nodes1[1]>>1) == (Nodes2[1]>>1))
        {
            continue;
        }



        if(asg_arc_a(g, Nodes1[1])[0].el == 0 || asg_arc_a(g, Nodes2[1])[0].el == 0)
        {
            continue;
        }

        uint32_t convex1, convex2, f1, f2;
        long long l1=-1, l2=-1;  // suppress compiler warning; l1 and l2 are always inited when used (by checking f1 and f2).
        todel = 0;

        if((Nodes1[0]^1) == (startNode^1) || (Nodes1[0]^1) == endNode)
        {
            continue;
        }
        if((Nodes2[1]^1) == (startNode^1) || (Nodes2[1]^1) == endNode)
        {
            continue;
        }

        f1 = detect_bubble_end_with_bubbles(g, Nodes1[0]^1, Nodes2[1]^1, &convex1, &l1, NULL);


        if((Nodes2[0]^1) == (startNode^1) || (Nodes2[0]^1) == endNode)
        {
            continue;
        }
        if((Nodes1[1]^1) == (startNode^1) || (Nodes1[1]^1) == endNode)
        {
            continue;
        }

        f2 = detect_bubble_end_with_bubbles(g, Nodes2[0]^1, Nodes1[1]^1, &convex2, &l2, NULL);


        if(f1 && f2)
        {
            if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
            {
                if(l1 <= min_thres || l2 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }
        else if(f1)
        {
           if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
            {
                if(l1 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }
        else if(f2)
        {
            if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
            {
                if(l2 <= min_thres)
                {
                    continue;
                }

                todel = 1;
            }
        }

        if(todel == 0)
        {
            if(!f1)
            {
                Ns_first[0] = Nodes1[0]^1; Ns_first[1] = Nodes2[1]^1;
                f1 = asg_bub_end_finder_with_del_advance(g, Ns_first, 2, max_dist, 
                bub, 0, (u_int32_t)-1, &convex1);
                l1 = min_thres + 10;
            }

            if(!f2)
            {
                Ns_second[0] = Nodes2[0]^1; Ns_second[1] = Nodes1[1]^1;
                f2 = asg_bub_end_finder_with_del_advance(g, Ns_second, 2, max_dist, 
                bub, 0, (u_int32_t)-1, &convex2);
                l2 = min_thres + 10;
            }

            if(f1 && f2)
            {
                if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)) && 
                ((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                {
                    if(l1 <= min_thres || l2 <= min_thres)
                    {
                        continue;
                    } 

                    todel = 1;
                } 
            }
            else if(f1)
            {
                if(((convex1>>1) == (startNode>>1) || (convex1>>1) == (endNode>>1)))
                {
                    if(l1 <= min_thres)
                    {
                        continue;
                    }

                    todel = 1;
                }   
            }
            else if(f2)
            {
                if(((convex2>>1) == (startNode>>1) || (convex2>>1) == (endNode>>1)))
                {
                    if(l2 <= min_thres)
                    {
                        continue;
                    }

                    todel = 1;
                }   
            }
            
        }

        
        if(todel)
        {
            if(av[0].el == 0 || aw[0].el == 0)
            {
                av[0].del = 1;
                asg_arc_del(g, av[0].v^1, av[0].ul>>32^1, 1);

                aw[0].del = 1;
                asg_arc_del(g, aw[0].v^1, aw[0].ul>>32^1, 1);
                
                n_reduced++;
            }
        }
    }

    return n_reduced;
}

int asg_arc_del_triangular_advance(asg_t *g, long long max_dist)
{
    double startTime = Get_T();
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, n_reduced_a = 0;

    if (!g->is_symm) asg_symm(g);

    buf_t b;  // bubble popping struct
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    buf_t bub;
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }


        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_triangular_exact(g, b.b.a, b.b.n, v, b.S.a[0], max_dist, &bub);
            n_reduced_a += test_triangular_addition_exact(g, b.b.a, b.b.n, v, b.S.a[0],max_dist, &bub);
        }

    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);

    if (n_reduced + n_reduced_a) {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d/%d triangular/triangular_a overlaps\n", 
        __func__, n_reduced, n_reduced_a);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced + n_reduced_a;
}






int check_if_cross(asg_t *g, uint32_t v)
{
    uint32_t N_list[5] = {0};
    if (g->seq[v>>1].del) return 0;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    if(nv != 2) return 0;
    if(asg_is_single_edge(g, av[0].v, v>>1) != 2 || asg_is_single_edge(g, av[1].v, v>>1) != 2)
    {
        return 0;
    }
    if(av[0].v == av[1].v)
    {
        return 0;
    }
    N_list[0] = v;
    N_list[1] = av[0].v^1;
    N_list[2] = av[1].v^1;
    if(asg_arc_n(g, N_list[0]) != 2 || 
       asg_arc_n(g, N_list[1]) != 2 ||
       asg_arc_n(g, N_list[2]) != 2 )
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[1])[0].v == asg_arc_a(g, N_list[1])[1].v)
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[2])[0].v == asg_arc_a(g, N_list[2])[1].v)
    {
        return 0;
    }

    if(asg_arc_a(g, N_list[1])[0].v == (N_list[0]^1))
    {
        N_list[3] = asg_arc_a(g, N_list[1])[1].v^1;
    }
    else if(asg_arc_a(g, N_list[1])[1].v == (N_list[0]^1))
    {
        N_list[3] = asg_arc_a(g, N_list[1])[0].v^1;
    }

    if(asg_arc_a(g, N_list[2])[0].v == (N_list[0]^1))
    {
        N_list[4] = asg_arc_a(g, N_list[2])[1].v^1;
    }
    else if(asg_arc_a(g, N_list[2])[1].v == (N_list[0]^1))
    {
        N_list[4] = asg_arc_a(g, N_list[2])[0].v^1;
    }

    if(N_list[3] != N_list[4])
    {
        return 0;
    }

    if(asg_arc_n(g, N_list[0]) != 2 || 
        asg_arc_n(g, N_list[1]) != 2 ||
        asg_arc_n(g, N_list[2]) != 2 ||
        asg_arc_n(g, N_list[3]) != 2)
    {
        return 0;
    }


    uint32_t convex1, convex2, f1, f2;
    long long l1, l2;
    l1 = l2 = 0;
    int todel = 0;
    
    f1 = detect_bubble_end_with_bubbles(g, N_list[0]^1, N_list[3]^1, &convex1, &l1, NULL);
    f2 = detect_bubble_end_with_bubbles(g, N_list[1]^1, N_list[2]^1, &convex2, &l2, NULL);

    if(f1 && f2)
    {
        if(l1 > min_thres && l2 > min_thres)
        {
            todel = 1;
        }     
    }
    else if(f1)
    {
        if(l1 > min_thres)
        {
            todel = 1;
        }        
    }
    else if(f2)
    {
        if(l2 > min_thres)
        {
            todel = 1;
        }  
    }
    

    return todel;
}




typedef struct {
	int threadID;
    int thread_num;
    int check_cross;
    asg_t *g;
} para_for_simple_bub;

void* asg_arc_identify_simple_bubbles_pthread(void* arg)
{
    int thr_ID = ((para_for_simple_bub*)arg)->threadID;
    int thr_num = ((para_for_simple_bub*)arg)->thread_num;
    asg_t *g = ((para_for_simple_bub*)arg)->g;
    int check_cross = ((para_for_simple_bub*)arg)->check_cross;
    ///the reason is that each read has two direction (query->target, target->query)    

	uint32_t v, w, n_vtx = g->n_seq * 2;
    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    long long l, i;
    ///for (v = 0; v < n_vtx; ++v)
    for (v = thr_ID; v < n_vtx; v = v + thr_num) 
    {
        if (g->seq[v>>1].del) continue;

        b.b.n = 0;  // b.b is an array of visited vertices (uint32_t)

        if(g->seq_vis[v] != 1)  // note: init was by calloc, to 0.
        {
            ///if(detect_bubble_with_bubbles(g, v, &w, &l, &b, (uint32_t)-1))
            if(detect_bubble_with_bubbles(g, v, &w, &l, &b, SMALL_BUBBLE_SIZE))
            {
                for (i = 0; i < (long long)b.b.n; i++)
                {
                    if(b.b.a[i] != (v>>1) && b.b.a[i] != (w>>1))
                    {
                        g->seq_vis[b.b.a[i]<<1] = 1;
                        g->seq_vis[(b.b.a[i]<<1) + 1] = 1;
                    }
                }
                g->seq_vis[v] = 1;
                g->seq_vis[w^1] = 1;
            }
        }

        if(check_cross == 1 && check_if_cross(g, v))
        {
            g->seq_vis[v] = 2;
        }
    }
    free(b.b.a);

    free(arg);

    return NULL;
}

int asg_arc_identify_simple_bubbles_multi(asg_t *g, int check_cross)
{
    double startTime = Get_T();
    memset(g->seq_vis, 0, g->n_seq*2*sizeof(uint8_t));

    pthread_t *_r_threads;

	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*asm_opt.thread_num);

    int i = 0;

	for (i = 0; i < asm_opt.thread_num; i++)
	{
        para_for_simple_bub* arg = (para_for_simple_bub*)malloc(sizeof(*arg));
        arg->g = g;
        arg->thread_num = asm_opt.thread_num;
        arg->threadID = i;
        arg->check_cross = check_cross;

        pthread_create(_r_threads + i, NULL, asg_arc_identify_simple_bubbles_pthread, (void*)arg);
	}
    

    for (i = 0; i<asm_opt.thread_num; i++)
		pthread_join(_r_threads[i], NULL);


    free(_r_threads);

    uint32_t v, n_vtx = g->n_seq * 2;
    long long nodes, bub_nodes, cross_nodes;
    bub_nodes = nodes = cross_nodes = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        nodes++;
        if(g->seq_vis[v] == 1) bub_nodes++;
        if(g->seq_vis[v] == 2) cross_nodes++;
    }


    fprintf(stderr, "[M::%s] dropped total of %lld, used %.1fs\n", 
                    __func__, bub_nodes+cross_nodes, Get_T()-startTime);
    return bub_nodes+cross_nodes;
}



int check_small_bubble(asg_t *g, uint32_t begNode, uint32_t v, uint32_t w, 
long long* vLen, long long* wLen, uint32_t* endNode)
{
    // my note: does not tolerate tips
    uint32_t nv = asg_arc_n(g, v);
    uint32_t nw = asg_arc_n(g, w);

    asg_arc_t *av = asg_arc_a(g, v);
    asg_arc_t *aw = asg_arc_a(g, w);
    if(nv != 1 || nw != 1)
    {
        return 0;
    }

    ///first node
    ///nv must be 1
    if(asg_is_single_edge(g, av[0].v, v>>1) == 2)  // the next node has 2 predecessors, potentially a bubble sink
    {
        uint32_t vv;
        vv = av[0].v^1;

        if(
            asg_is_single_edge(g, asg_arc_a(g, vv)[0].v, vv>>1) == 1
            &&
            asg_is_single_edge(g, asg_arc_a(g, vv)[1].v, vv>>1) == 1
            )
        {
            ///walk along first path
            //  note: pLen1 will be -1 if no single path is available between `asg_arc_a(g, vv)[0].v` and `begNode>>1`
            long long pLen1;
            pLen1 = single_edge_length(g, asg_arc_a(g, vv)[0].v, begNode>>1, 1000);
            long long pLen2;
            pLen2 = single_edge_length(g, asg_arc_a(g, vv)[1].v, begNode>>1, 1000);

            if(pLen1 >= 0 && pLen2 >= 0)
            {
                if(((asg_arc_a(g, vv)[0].v) == (v^1)) && pLen1 == 1)
                {
                    (*vLen) = pLen1;
                    (*wLen) = pLen2;
                }
                else if(((asg_arc_a(g, vv)[1].v) == (v^1)) && pLen2 == 1)
                {
                    (*vLen) = pLen2;
                    (*wLen) = pLen1;
                }
                else
                {
                    fprintf(stderr, "[E::%s]error\n", __func__);
                }
                ///(*endNode) = vv>>1;
                (*endNode) = vv;
                return 1;
            }

        }
    }


    ///second node
    ///nw must be 1
    if(asg_is_single_edge(g, aw[0].v, w>>1) == 2)
    {
        uint32_t ww;
        ww = aw[0].v^1;

        if(
            asg_is_single_edge(g, asg_arc_a(g, ww)[0].v, ww>>1) == 1
            &&
            asg_is_single_edge(g, asg_arc_a(g, ww)[1].v, ww>>1) == 1
            )
        {
            ///walk along first path
            long long pLen1;
            pLen1 = single_edge_length(g, asg_arc_a(g, ww)[0].v, begNode>>1, 1000);

            ///walk along first path
            long long pLen2;
            pLen2 = single_edge_length(g, asg_arc_a(g, ww)[1].v, begNode>>1, 1000);

            if(pLen1 >= 0 && pLen2 >= 0)
            {
                if(((asg_arc_a(g, ww)[0].v) == (w^1)) && pLen1 == 1)
                {
                    (*wLen) = pLen1;
                    (*vLen) = pLen2;
                }
                else if(((asg_arc_a(g, ww)[1].v) == (w^1)) && pLen2 == 1)
                {
                    (*wLen) = pLen2;
                    (*vLen) = pLen1;
                }
                else
                {
                    fprintf(stderr, "[E::%s]error\n", __func__);
                }

                //(*endNode) = ww>>1;
                (*endNode) = ww;

                return 1;
            }
       }
    }

    return 0;

}


int test_single_node_bubble(asg_t *g, uint32_t* nodes, uint32_t length, 
uint32_t startNode, uint32_t endNode)
{
    
    uint32_t i, v, w;
    uint32_t vEnd;
    int flag0, flag1;
    int n_reduced = 0;
    long long Len[2], longLen;
    long long longLen_thres = 4;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }

        
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);
        if(flag0 != 1 || flag1 != 1)
        {
            continue;
        }

        if(check_small_bubble(g, v, av[0].v, av[1].v, &(Len[0]), &(Len[1]), &vEnd))
        {

            if(if_node_exist(nodes, length, vEnd>>1) && ((vEnd>>1) != (endNode>>1)))
           {
               if(Len[0] == 1 && Len[1] != 1)  
               {
                   w = av[0].v;
                   longLen = Len[1];
                    
                  /****************************may have bugs********************************/
                  if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                   {
                    //    fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                    //    w>>1, startNode>>1, endNode>>1);
                       asg_seq_del(g, w>>1);
                       n_reduced++;
                   }///up to here w is exactly overlapped in both directions
                   else if(longLen >= longLen_thres)
                   {
                       if(av[0].el == 1 && av[1].el == 1
                          &&
                          asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                       {
                           asg_seq_del(g, w>>1);
                           n_reduced++;
                       }
                   }
                   
                  /****************************may have bugs********************************/
                   
               }
               else if(Len[0] != 1 && Len[1] == 1)
               {
                   w = av[1].v;
                   longLen = Len[0];
                   
                  /****************************may have bugs********************************/
                  if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                   {
                    //    fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                    //    w>>1, startNode>>1, endNode>>1);
                       asg_seq_del(g, w>>1);
                       n_reduced++;
                   }///up to here w is exactly overlapped in both directions
                   else if(longLen >= longLen_thres)
                   {
                       if(av[0].el == 1 && av[1].el == 1
                          &&
                          asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                       {
                           asg_seq_del(g, w>>1);
                           n_reduced++;
                       }
                   }
                  /****************************may have bugs********************************/
               }
               else if(Len[0] == 1 && Len[1] == 1)
               {
                   w = av[0].v;
                   flag0 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                   w = av[1].v;
                   flag1 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                    ///>=2 means this is an exact overlap
                   if(flag0 < 2 && flag1 >= 2)
                   {
                       w = av[0].v;
                       
                       /****************************may have bugs********************************/
                        if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                        {
                            // fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                            // w>>1, startNode>>1, endNode>>1);
                            asg_seq_del(g, w>>1);
                            n_reduced++;
                        }
                        /****************************may have bugs********************************/
                   }

                   if(flag0 >= 2 && flag1 < 2)
                   {
                       w = av[1].v;
                      
                       /****************************may have bugs********************************/
                        if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                        {
                            // fprintf(stderr, "w>>1: %u, beg: %u, end: %u\n", 
                            // w>>1, startNode>>1, endNode>>1);
                            asg_seq_del(g, w>>1);
                            n_reduced++;
                        }
                        /****************************may have bugs********************************/
                   }
               }
               else
               {
                   fprintf(stderr, "[E::%s]error\n", __func__);
               }
               
           }
        }

    }

    return n_reduced;
}


int test_single_node_bubble_directly(asg_t *g, uint32_t v, long long longLen_thres, ma_hit_t_alloc* sources)
{
    uint32_t w, vEnd;
    int flag0, flag1;
    int n_reduced = 0;
    long long Len[2], longLen;
    ///long long longLen_thres = 4;

    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    if(nv != 2)
    {
        return 0;
    }

    
    flag0 = asg_is_single_edge(g, av[0].v, v>>1);  // function name does not say what it does; count of predecessor of the given vertex
    flag1 = asg_is_single_edge(g, av[1].v, v>>1);  // function name does not say what it does; count of predecessor of the given vertex
    if(flag0 != 1 || flag1 != 1)  // not a simplest start of a bubble
    {
        return 0;
    }

    if(check_small_bubble(g, v, av[0].v, av[1].v, &(Len[0]), &(Len[1]), &vEnd))
    {
        // ~remove the shroter side~
        if(Len[0] == 1 && Len[1] != 1)
        {
            w = av[0].v;
            longLen = Len[1];
            
            /****************************may have bugs********************************/
            ///if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
            if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0 || sources[w>>1].is_abnormal == 1)
            {
                asg_seq_del(g, w>>1);
                n_reduced++;
            }///up to here w is exactly overlapped in both directions
            else if(longLen >= longLen_thres)
            {
                if(av[0].el == 1 && av[1].el == 1
                    &&
                    asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                {
                    asg_seq_del(g, w>>1);
                    n_reduced++;
                }
            }
            
            /****************************may have bugs********************************/
            
        }
        else if(Len[0] != 1 && Len[1] == 1)
        {
            w = av[1].v;
            longLen = Len[0];
            
            /****************************may have bugs********************************/
            ///if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
            if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0 || sources[w>>1].is_abnormal == 1)
            {
                asg_seq_del(g, w>>1);
                n_reduced++;
            }///up to here w is exactly overlapped in both directions
            else if(longLen >= longLen_thres)
            {
                if(av[0].el == 1 && av[1].el == 1
                    &&
                    asg_arc_a(g, vEnd)[0].el == 1 && asg_arc_a(g, vEnd)[1].el == 1)
                {
                    asg_seq_del(g, w>>1);
                    n_reduced++;
                }
            }
            /****************************may have bugs********************************/
        }
        else if(Len[0] == 1 && Len[1] == 1)  // both sides are equally long: 1 arc away
        {
            flag0 = sources[av[0].v>>1].is_abnormal;
            flag1 = sources[av[1].v>>1].is_abnormal;

            if(flag0 == 1 && flag1 == 0)
            {
                asg_seq_del(g, av[0].v>>1);
                n_reduced++;
            }
            else if(flag0 == 0 && flag1 == 1)
            {
                asg_seq_del(g, av[1].v>>1);
                n_reduced++;
            }
            else
            {
                w = av[0].v;
                flag0 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                w = av[1].v;
                flag1 = asg_arc_a(g, w)[0].el + asg_arc_a(g, w^1)[0].el;
                ///>=2 means this is an exact overlap
                if(flag0 < 2 && flag1 >= 2)
                {
                    w = av[0].v;
                    
                    /****************************may have bugs********************************/
                    if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                    {
                        asg_seq_del(g, w>>1);
                        n_reduced++;
                    }
                    /****************************may have bugs********************************/
                }

                if(flag0 >= 2 && flag1 < 2)
                {
                    w = av[1].v;
                    
                    /****************************may have bugs********************************/
                    if(asg_arc_a(g, w)[0].el == 0 || asg_arc_a(g, w^1)[0].el == 0)
                    {
                        asg_seq_del(g, w>>1);
                        n_reduced++;
                    }
                    /****************************may have bugs********************************/
                }
            }
        }
        
    }
    return n_reduced;
}


int asg_arc_del_single_node_bubble(asg_t *g, long long max_dist)
{
    ///the reason is that each read has two direction (query->target, target->query)
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_single_node_bubble(g, b.b.a, b.b.n, v, b.S.a[0]);
        }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d short bubbles\n\n", __func__, n_reduced);

    return n_reduced;
}

int asg_arc_del_single_node_directly(asg_t *g, long long longLen_thres, ma_hit_t_alloc* sources)
{  // pop small bubble
   // if one side is 1-arc long and the other side appears longer, the short side's node is removed
   // if both 1-arc long, it depends on alignment
   // (doesn't process larger/more complex bubbles.)
    double startTime = Get_T();
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv != 2)
        {
            continue;
        }

        n_reduced += test_single_node_bubble_directly(g, v, longLen_thres, sources);
    }


    if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    return n_reduced;
}



int test_cross(asg_t *g, uint32_t* nodes, uint32_t length, 
               uint32_t startNode, uint32_t endNode)
{
    uint32_t a1, a2;
    uint32_t N_list[5] = {0};
    uint32_t i, v;
    int flag0, flag1;
    int n_reduced = 0;
    for (i = 0; i < length; ++i)
    {
        v = nodes[i];

        if((v>>1) == (startNode>>1) || (v>>1) == (endNode>>1))
        {
            continue;
        }

        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if(nv != 2)
        {
            continue;
        }
        if(av[0].v == av[1].v)
        {
            continue;
        }
        flag0 = asg_is_single_edge(g, av[0].v, v>>1);
        flag1 = asg_is_single_edge(g, av[1].v, v>>1);

        if(flag0 != 2 || flag1 != 2)
        {
            continue;
        }
        

        N_list[0] = v;
        N_list[1] = av[0].v^1;
        N_list[2] = av[1].v^1;

        if(asg_arc_n(g, N_list[0]) != 2 || 
           asg_arc_n(g, N_list[1]) != 2 ||
           asg_arc_n(g, N_list[2]) != 2 )
        {
            continue;
        }

        if(asg_arc_a(g, N_list[1])[0].v == asg_arc_a(g, N_list[1])[1].v)
        {
            continue;
        }

        if(asg_arc_a(g, N_list[2])[0].v == asg_arc_a(g, N_list[2])[1].v)
        {
            continue;
        }


        if(asg_arc_a(g, N_list[1])[0].v == (N_list[0]^1))
        {
            N_list[3] = asg_arc_a(g, N_list[1])[1].v^1;
        }
        else if(asg_arc_a(g, N_list[1])[1].v == (N_list[0]^1))
        {
            N_list[3] = asg_arc_a(g, N_list[1])[0].v^1;
        }
        else
        {
            fprintf(stderr, "[E::%s]ERROR at %s:%d\n", __func__, __FILE__, __LINE__);
        }
        
        if(asg_arc_a(g, N_list[2])[0].v == (N_list[0]^1))
        {
            N_list[4] = asg_arc_a(g, N_list[2])[1].v^1;
        }
        else if(asg_arc_a(g, N_list[2])[1].v == (N_list[0]^1))
        {
            N_list[4] = asg_arc_a(g, N_list[2])[0].v^1;
        }
        else
        {
            fprintf(stderr, "[E::%s]ERROR at %s:%d\n", __func__, __FILE__, __LINE__);
        }

        if(N_list[3] != N_list[4])
        {
            continue;
        }

        if(asg_arc_n(g, N_list[0]) != 2 || 
           asg_arc_n(g, N_list[1]) != 2 ||
           asg_arc_n(g, N_list[2]) != 2 ||
           asg_arc_n(g, N_list[3]) != 2)
        {
            continue;
        }
        /**
        N_list[3]        N_list[0]

        N_list[2]        N_list[1]
        **/
       if(asg_arc_a(g, N_list[0])[0].el == asg_arc_a(g, N_list[0])[1].el)
       {
           continue;
       }

       if(asg_arc_a(g, N_list[0])[0].el == 1)
       {
           //a1 = asg_arc_a(g, N_list[0])[0].v >> 1;
           a1 = 0;
       }
       else
       {
           ///a1 = asg_arc_a(g, N_list[0])[1].v >> 1;
           a1 = 1;
       }





       
       if(asg_arc_a(g, N_list[3])[0].el == asg_arc_a(g, N_list[3])[1].el)
       {
           continue;
       }

       if(asg_arc_a(g, N_list[3])[0].el == 1)
       {
           //a2 = asg_arc_a(g, N_list[3])[0].v >> 1;
           a2 = 0;
       }
       else
       {
           //a2 = asg_arc_a(g, N_list[3])[1].v >> 1;
           a2 = 1;
       }

        if(
            (asg_arc_a(g, N_list[0])[a1].v >> 1)
            !=
            (asg_arc_a(g, N_list[3])[a2].v >> 1)
           )
        {
            if(((N_list[0]>>1) != (endNode>>1)) &&
               ((N_list[1]>>1) != (endNode>>1)) && 
               ((N_list[2]>>1) != (endNode>>1)) &&
               ((N_list[3]>>1) != (endNode>>1)))
            {
                asg_arc_a(g, N_list[0])[a1].del = 1;
                asg_arc_del(g, asg_arc_a(g, N_list[0])[a1].v^1, 
                asg_arc_a(g, N_list[0])[a1].ul>>32^1, 1);


                asg_arc_a(g, N_list[3])[a2].del = 1;
                asg_arc_del(g, asg_arc_a(g, N_list[3])[a2].v^1, 
                asg_arc_a(g, N_list[3])[a2].ul>>32^1, 1);
                n_reduced++;
            }
        }
    }

    return n_reduced;
}

int asg_arc_del_cross_bubble(asg_t *g, long long max_dist)
{
    double startTime = Get_T();
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }

        ///if this is a bubble
        if(asg_bub_finder_with_del_advance(g, v, max_dist, &b) == 1)
        {
            n_reduced += test_cross(g, b.b.a, b.b.n, v, b.S.a[0]);
        }
    
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);

    if (n_reduced) {        
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d cross\n", __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}



// transitive reduction; see Myers, 2005
int asg_arc_del_trans(asg_t *g, int fuzz)
{
    double startTime = Get_T();

	uint8_t *mark;
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
	///at first, all nodes should be set to vacant
	mark = (uint8_t*)calloc(n_vtx, 1);

	/**v is the id+direction of a node, 
     * the high 31-bit is the id, 
     * and the lowest 1-bit is the direction
     * (0 means query-to-target, 1 means target-to-query)**/
	for (v = 0; v < n_vtx; ++v) {
		uint32_t L, i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);

        // current node is removed, or the vertex has no target
		if (nv == 0) continue;
		if (g->seq[v>>1].del) 
        {
			for (i = 0; i < nv; ++i) av[i].del = 1, ++n_reduced;
			continue;
		}

        /**
	********************************query-to-target overlap****************************
	case 1: u = 0, rev = 0                           in the view of target: direction is 1 
	query: CCCCCCCCTAATTAAAAT                        target: TAATTAAAATGGGGGG (use ex-target as query)
	               ||||||||||         <--->                  ||||||||||
	       target: TAATTAAAATGGGGGG           query: CCCCCCCCTAATTAAAAT (use ex-query as target)

	case 2: u = 0, rev = 1                           in the view of target: direction is 0
	query: CCCCCCCCTAATTAAAAT					     target: CCCCCCATTTTAATTA  (use ex-target as query)
                   ||||||||||        <--->                         ||||||||||
	       target: TAATTAAAATGGGGGG                         query: ATTTTAATTAGGGGGGGG  (use ex-query as target)
	********************************query-to-target overlap****************************

	********************************target-to-query overlap****************************
	case 3: u = 1, rev = 0                           in the view of target: direction is 0
			 query: AAATAATATCCCCCCGCG                target: GGGCCGGCAAATAATAT (use ex-target as query)
					|||||||||          <--->                          |||||||||
	target: GGGCCGGCAAATAATAT                                  query: AAATAATATCCCCCCGCG (use ex-query as target)

	case 4: u = 1, rev = 1                          in the view of target: direction is 1
	         query: AAATAATATCCCCCCGCG                        target: ATATTATTTGCCGGCCC (use ex-target as query)
                    |||||||||             <--->                       |||||||||
	target: GGGCCGGCAAATAATAT                          query: CGCGGGGGATATTATTT (use ex-query as target)
	********************************target-to-query overlap****************************

	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
						                (in the view of query)
	p->v : |___________31___________|__________1___________|
				        tns             reverse direction of overlap
						              (in the view of target)
	p->ol: overlap length
    **/


		//all outnode of v should be set to "not reduce"
		for (i = 0; i < nv; ++i) mark[av[i].v] = 1;

		///length of node (not overlap length)
		///av[nv-1] is longest out-dege
		/**
		 * v---------------
		 *   w1---------------
		 *      w2--------------
		 *         w3--------------
		 *            w4--------------
		 *               w5-------------
		 * for v, the longest out-edge is v->w5
		 **/
		L = asg_arc_len(av[nv-1]) + fuzz;


		for (i = 0; i < nv; ++i) {
			//w is an out-node of v
			uint32_t w = av[i].v;
			
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			///if w has already been reduced
			if (mark[av[i].v] != 1) continue;

			for (j = 0; j < nw && asg_arc_len(aw[j]) + asg_arc_len(av[i]) <= L; ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#if 0
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || asg_arc_len(aw[j]) < fuzz); ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#endif
        //remove edges
		for (i = 0; i < nv; ++i) {
			if (mark[av[i].v] == 2) av[i].del = 1, ++n_reduced;
			mark[av[i].v] = 0;
		}
	}
	free(mark);

	if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] reduced %d arcs, used %.1fs\n", 
                    __func__, n_reduced, Get_T()-startTime);
	return n_reduced;
}





///max_ext is 4
int asg_cut_tip(asg_t *g, int max_ext)
{
    double startTime = Get_T();

	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    
	for (v = 0; v < n_vtx; ++v) {
		if (g->seq[v>>1].del) continue;
		///check if the another direction of v has no overlaps
		///if the self direction of v has no overlaps, we don't have the overlaps of them
		///here is check if the reverse direction of v 
		/**
		 the following first line is to find (means v is a node has no prefix):
		        (v)--->()---->()---->()----->....
         another case is:
            ......()---->()---->()----->()------>(v)
         this case can be found by (v^1), so we don't need to process this case here
		**/
		if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
        /**
         the following second line is: 
		 (v)--->()---->()---->()----->()
		        |--------max_ext-------|
        **/
       ///that means here is a long tip, which is longer than max_ext
		if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig

		/**
		 * so combining the last two lines, they are designed to reomve(n(0), n(1), n(2)):
		 *                  		    ----->n(4)
		 *                  		   |
		 * n(0)--->n(1)---->n(2)---->n(3)
		 *                             |
		 *                              ----->n(5)
		**/
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
    
    fprintf(stderr, "[M::%s] cut %d tips, used %.1f s\n", 
                    __func__, cnt, Get_T()-startTime);
	return cnt;
}


///max_ext is 4
int asg_cut_tip_primary(asg_t *g, ma_ug_t *ug, int max_ext)
{
    double startTime = Get_T();

	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0, tipEvaluateLen;
    
	for (v = 0; v < n_vtx; ++v) {
		//if this seq has been deleted
		if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
		///(v^1) is end node of a unitig
		if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
        /**
         the following second line is: 
		 (v)--->()---->()---->()----->()
		        |--------max_ext-------|
        **/
       ///that means here is a long tip, which is longer than max_ext
		if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig

		/**
		 * so combining the last two lines, they are designed to reomve(n(0), n(1), n(2)):
		 *                  		    ----->n(4)
		 *                  		   |
		 * n(0)--->n(1)---->n(2)---->n(3)
		 *                             |
		 *                              ----->n(5)
		**/

        if(ug!=NULL)
        {
            tipEvaluateLen = 0;
            for (i = 0; i < a.n; ++i)
            {
                tipEvaluateLen += EvaluateLen(ug->u, ((uint32_t)a.a[i]>>1));
            }
            if(tipEvaluateLen > (uint32_t)max_ext) continue;
        }
    
        for (i = 0; i < a.n; ++i)
        {
            g->seq[((uint32_t)a.a[i]>>1)].c = ALTER_LABLE;
        }

        for (i = 0; i < a.n; ++i)
        {
            asg_seq_drop(g, (uint32_t)a.a[i]>>1);
        }

		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return cnt;
}





// delete short arcs
///for best graph?
int asg_arc_del_short(asg_t *g, float drop_ratio)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
	for (v = 0; v < n_vtx; ++v) {
		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
	}
	if (n_short) {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	return n_short;
}


inline int check_weak_ma_hit(ma_hit_t_alloc* aim_paf, ma_hit_t_alloc* reverse_paf_list, 
long long weakID, uint32_t w_qs, uint32_t w_qe)
{
    long long i = 0;
    long long strongID, index;
    for (i = 0; i < aim_paf->length; i++)
    {
        ///if this is a strong overlap
        if (
        aim_paf->buffer[i].del == 0
        &&
        aim_paf->buffer[i].ml == 1
        && 
        Get_qs(aim_paf->buffer[i]) <= w_qs
        && 
        Get_qe(aim_paf->buffer[i]) >= w_qe)
        {
            strongID = Get_tn(aim_paf->buffer[i]);
            index = get_specific_overlap(&(reverse_paf_list[strongID]), strongID, weakID);
            if(index != -1)
            {
                return 0;
            }
        }
    }

    return 1;
}


inline int check_weak_ma_hit_reverse(ma_hit_t_alloc* r_paf, ma_hit_t_alloc* r_paf_source, 
long long weakID)
{
    long long i = 0;
    long long strongID, index;
    ///all overlaps coming from another haplotye are strong
    for (i = 0; i < r_paf->length; i++)
    {
        strongID = Get_tn(r_paf->buffer[i]);
        index = get_specific_overlap
            (&(r_paf_source[strongID]), strongID, weakID);
        ///must be a strong overlap
        if(index != -1 && r_paf_source[strongID].buffer[index].ml == 1)
        {
            return 0;
        }
    }

    return 1;
}


inline int check_weak_ma_hit_debug(ma_hit_t_alloc* aim_paf, ma_hit_t_alloc* reverse_paf_list, 
long long weakID)
{
    long long i = 0;
    long long strongID, index;
    for (i = 0; i < aim_paf->length; i++)
    {
        ///if this is a strong overlap
        if (aim_paf->buffer[i].ml == 1)
        {
            strongID = Get_tn(aim_paf->buffer[i]);
            index = get_specific_overlap(&(reverse_paf_list[strongID]), strongID, weakID);
            if(index != -1)
            {
                return strongID;
            }
        }


    }

    return 0;
}



// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_unclean(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources)
{
    double startTime = Get_T();

	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    uint32_t last_e;
	for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;

		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		//av[0] has the most overlap length
		///remove short overlaps
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		///av has been sorted by overlap length
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i) {}
        last_e = i + 1;

		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
        
        
        if(nv >= 2 && av[1].del == 1)
        {
            ///second longest
            av[1].del = 0;
            --n_short;
            last_e++;
        }
        
	}
	///if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

	return n_short;
}









/*****************************read graph*****************************/

uint32_t detect_single_path_with_dels(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b)
{
    // FUNC
    //     walk the single path start from begNode
    //     (theres no limit of the single path length)
    // SIDE EFFECT
    //     the node before branching / is the end in *endNode
    //     path length (in number of nodes) in *Len
    //     walked nodes in *b
    // RETURN
    //     (see macros)
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;

    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            return END_TIPS;
        }

        if(kv == 2)
        {
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

uint32_t detect_single_path_with_dels_n_stops(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* Len, long long* max_stop_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw, n_stops = 0, flag = LONG_TIPS;
    (*Len) = 0;
    (*max_stop_Len) = 0;
    long long preLen = 0, currentLen;

    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            flag = MUL_OUTPUT;
            break;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;

        if(kw >= 2)
        {
            n_stops++;
            currentLen = (*Len) - preLen;
            preLen = (*Len);
            if(currentLen > (*max_stop_Len))
            {
                (*max_stop_Len) = currentLen;
            }
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            (*Len)++;
            if(b) kv_push(uint32_t, b->b, v>>1);
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        } 
    }


    currentLen = (*Len) - preLen;
    preLen = (*Len);
    if(currentLen > (*max_stop_Len))
    {
        (*max_stop_Len) = currentLen;
    }
    return flag;   
}

uint32_t detect_single_path_with_dels_contigLen(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* baseLen, buf_t* b)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k;
    (*baseLen) = 0;


    while (1)
    {
        ///(*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            return END_TIPS;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return MUL_OUTPUT;
        }

        ///kv must be 1
        for (k = 0; k < asg_arc_n(g, v); k++)
        {
            if(!asg_arc_a(g, v)[k].del)
            {
                w = asg_arc_a(g, v)[k].v;
                (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                break;
            }
        }

        ///up to here, kv=1
        ///kw must >= 1
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

uint32_t detect_single_path_with_dels_contigLen_complex(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, long long* max_stop_base_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k, n_stops = 0, flag = LONG_TIPS;
    (*baseLen) = 0;
    (*max_stop_base_Len) = 0;
    long long preBaseLen = 0, currentBaseLen;


    while (1)
    {
        ///(*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = MUL_OUTPUT;
            break;
        }

        ///kv must be 1
        for (k = 0; k < asg_arc_n(g, v); k++)
        {
            if(!asg_arc_a(g, v)[k].del)
            {
                w = asg_arc_a(g, v)[k].v;
                (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                break;
            }
        }

        ///up to here, kv=1
        ///kw must >= 1
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;

        if(kw >= 2)
        {
            n_stops++;
            currentBaseLen = (*baseLen) - preBaseLen;
            preBaseLen = (*baseLen);
            if(currentBaseLen > (*max_stop_base_Len))
            {
                (*max_stop_base_Len) = currentBaseLen;
            }
        }

        if(kw >= 2 && n_stops >= stops_threshold)
        {
            (*baseLen) += g->seq[v>>1].len;
            if(b) kv_push(uint32_t, b->b, v>>1);
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }


        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        } 
    }

    currentBaseLen = (*baseLen) - preBaseLen;
    preBaseLen = (*baseLen);
    if(currentBaseLen > (*max_stop_base_Len))
    {
        (*max_stop_base_Len) = currentBaseLen;
    }

    return flag;   
}

/*****************************read graph*****************************/




/*****************************untig graph*****************************/

uint32_t untig_detect_single_path_with_dels(asg_t *g, ma_ug_t *ug, uint32_t begNode, 
uint32_t* endNode, long long* Len, buf_t* b)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;


    while (1)
    {
        (*Len) = (*Len) + EvaluateLen((*ug).u, v>>1);
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            return END_TIPS;
        }

        if(kv == 2)
        {
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            return MUL_OUTPUT;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);

        if(kw == 2)
        {
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            return MUL_INPUT;
        }   

        v = w;

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
    }

    return LONG_TIPS;   
}

uint32_t untig_detect_single_path_with_dels_n_stops(asg_t *g, ma_ug_t *ug, uint32_t begNode, 
uint32_t* endNode, long long* Len, long long* max_stop_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw, n_stops = 0, flag = LONG_TIPS;
    (*Len) = 0;
    (*max_stop_Len) = 0;
    long long preLen = 0, currentLen;

    while (1)
    {
        (*Len) = (*Len) + EvaluateLen((*ug).u, v>>1);
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            flag = MUL_OUTPUT;
            break;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);

        ///just calculate the max_stop_Len
        if(kw >= 2)
        {
            n_stops++;
            currentLen = (*Len) - preLen;
            preLen = (*Len);
            if(currentLen > (*max_stop_Len))
            {
                (*max_stop_Len) = currentLen;
            }
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        v = w; 

        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        }
    }


    currentLen = (*Len) - preLen;
    preLen = (*Len);
    if(currentLen > (*max_stop_Len))
    {
        (*max_stop_Len) = currentLen;
    }
    return flag;   
}

uint32_t untig_detect_single_path_with_dels_contigLen(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, buf_t* b)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k;
    (*baseLen) = 0;


    while (1)
    {
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            return END_TIPS;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return MUL_OUTPUT;
        }


        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);
        
        if(kw == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            return MUL_INPUT;
        }


        for (k = 0; k < asg_arc_n(g, v); k++)
        {
            if(!asg_arc_a(g, v)[k].del)
            {
                w = asg_arc_a(g, v)[k].v;
                (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                break;
            }
        }


        v = w; 

        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        }
    }

    return LONG_TIPS;   
}

uint32_t untig_detect_single_path_with_dels_contigLen_complex(asg_t *g, uint32_t begNode, uint32_t* endNode, 
long long* baseLen, long long* max_stop_base_Len, buf_t* b, uint32_t stops_threshold)
{
    
    uint32_t v = begNode, w = 0;
    uint32_t kv, kw, k, n_stops = 0, flag = LONG_TIPS;
    (*baseLen) = 0;
    (*max_stop_base_Len) = 0;
    long long preBaseLen = 0, currentBaseLen;


    while (1)
    {
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        if(b) kv_push(uint32_t, b->b, v>>1);

        if(kv == 0)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = END_TIPS;
            break;
        }

        if(kv == 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = TWO_OUTPUT;
            break;
        }

        if(kv > 2)
        {
            (*baseLen) += g->seq[v>>1].len;
            flag = MUL_OUTPUT;
            break;
        }

        

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);
        kw = get_real_length(g, w^1, NULL);

        if(kw == 1)
        {
            ///kv must be 1
            for (k = 0; k < asg_arc_n(g, v); k++)
            {
                if(!asg_arc_a(g, v)[k].del)
                {
                    w = asg_arc_a(g, v)[k].v;
                    (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                    break;
                }
            }
        }

        
        ///just calculate the max_stop_Len
        if(kw >= 2)
        {
            n_stops++;
            if(n_stops >= stops_threshold)
            {
                (*baseLen) += g->seq[v>>1].len;
            }
            else
            {
                for (k = 0; k < asg_arc_n(g, v); k++)
                {
                    if(!asg_arc_a(g, v)[k].del)
                    {
                        w = asg_arc_a(g, v)[k].v;
                        (*baseLen) += ((uint32_t)(asg_arc_a(g, v)[k].ul));
                        break;
                    }
                }
            }
            
            currentBaseLen = (*baseLen) - preBaseLen;
            preBaseLen = (*baseLen);
            if(currentBaseLen > (*max_stop_base_Len))
            {
                (*max_stop_base_Len) = currentBaseLen;
            }
        }

        if(kw >= 2 && n_stops >= stops_threshold)
        {
            if(kw == 2) flag = TWO_INPUT;
            if(kw > 2) flag = MUL_INPUT;
            break;
        }
        
        
        v = w; 

        if((v>>1) == (begNode>>1))
        {
            flag = LOOP;
            break;
        }
    }

    currentBaseLen = (*baseLen) - preBaseLen;
    preBaseLen = (*baseLen);
    if(currentBaseLen > (*max_stop_base_Len))
    {
        (*max_stop_base_Len) = currentBaseLen;
    }

    return flag;   
}

/*****************************untig graph*****************************/







long long check_if_diploid(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, R_to_U* ruIndex)
{
    // RETURN
    //     -1 if either v1 or v2 found loop in their single paths, or the single paths are too short
    //     
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
    long long l1, l2;

    b_0.b.n = 0;
    b_1.b.n = 0;
    ///uint32_t flag1 = detect_single_path(g, v1, &convex1, &l1, &b_0);
    uint32_t flag1 = detect_single_path_with_dels(g, v1, &convex1, &l1, &b_0);
    
    ///uint32_t flag2 = detect_single_path(g, v2, &convex2, &l2, &b_1);
    uint32_t flag2 = detect_single_path_with_dels(g, v2, &convex2, &l2, &b_1);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
		free(b_0.b.a); free(b_1.b.a);
        return -1;
    }

    // adjust the stored path length
    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }
    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }

    // single path too short
    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
		free(b_0.b.a); free(b_1.b.a);
        return -1;
    }

    buf_t* b_min;  // nodes in the shorter walk
    buf_t* b_max;  // nodes in the longer walk
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }

    long long i, j, k;
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn, is_Unitig;
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            /****************************may have bugs********************************/
            ///if(g->seq[tn].del == 1) continue;
            if(g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || g->seq[tn].del == 1) continue;
            }
            /****************************may have bugs********************************/
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }


    free(b_0.b.a);
    free(b_1.b.a);

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count/min_count>0.3) return 1;
    return 0;

}


long long check_if_diploid_primary_complex(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, uint32_t stops_threshold,
int if_drop, R_to_U* ruIndex)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
    long long l1, l2, max_stop_Len;

    b_0.b.n = 0;
    b_1.b.n = 0;
    uint32_t flag1 = detect_single_path_with_dels_n_stops(g, v1, &convex1, &l1, 
    &max_stop_Len, &b_0, stops_threshold);
    
    uint32_t flag2 = detect_single_path_with_dels_n_stops(g, v2, &convex2, &l2, 
    &max_stop_Len, &b_1, stops_threshold);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
        return -1;
    }

    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }

    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }


    long long i, j, k;
    i = b_0.b.n; i--;
    j = b_1.b.n; j--;
    while (i>=0 && j>=0)
    {
        if(b_0.b.a[i] == b_1.b.a[j])
        {
            i--;
            j--;
        }
        else
        {
            break;
        }
        
    }
    b_0.b.n = i+1; l1 = i+1;
    b_1.b.n = j+1; l2 = j+1;


    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_min;
    buf_t* b_max;
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }

    
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn, is_Unitig;
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            /****************************may have bugs********************************/
            ///if(g->seq[tn].del == 1 || (if_drop == 1 && g->seq[tn].c == ALTER_LABLE)) continue;
            if(g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || g->seq[tn].del == 1) continue;
            }
            /****************************may have bugs********************************/
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }


    free(b_0.b.a);
    free(b_1.b.a);
    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count/min_count>0.3) return 1;
    return 0;

}


uint32_t get_long_tip_length_stops(asg_t *sg, ma_utg_v* u, uint32_t begNode, uint32_t* endNode, 
buf_t* b, uint32_t stops_threshold)
{
    uint32_t v = begNode, w, n_stops = 0;
    uint32_t kv, kw;
    uint32_t eLen = 0;
    (*endNode) = (uint32_t)-1;
    if(u->a[v>>1].circ) return eLen;
    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
        eLen += EvaluateLen((*u), v>>1);
        if(b) kv_push(uint32_t, b->b, v);
        if(kv!=1) return eLen;
        ///kw must be 1 here
        kw = get_real_length(sg, v, &w);
        kw = get_real_length(sg, w^1, NULL);
        ///if(get_real_length(sg, w^1, NULL)!=1) return eLen;
        if(kw >= 2)
        {
            n_stops++;
        }
        if(kw >= 2 && n_stops >= stops_threshold)
        {
            return eLen;
        }
        v = w;
        if(v == begNode) return eLen;
    }
}



uint32_t get_long_tip_length(asg_t *sg, ma_utg_v* u, uint32_t begNode, uint32_t* endNode, buf_t* b)
{
    uint32_t v = begNode, w;
    uint32_t kv;
    uint32_t eLen = 0;
    (*endNode) = (uint32_t)-1;
    if(u->a[v>>1].circ) return eLen;
    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        (*endNode) = v;
        eLen += EvaluateLen((*u), v>>1);
        if(b) kv_push(uint32_t, b->b, v);
        if(kv!=1) return eLen;
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
        if(get_real_length(sg, w^1, NULL)!=1) return eLen;
        v = w;
        if(v == begNode)
        {
            u->a[begNode>>1].circ = 1;
            return eLen;
        } 
    }
}

uint32_t if_long_tip_length(asg_t *sg, ma_ug_t *ug, uint32_t begNode, uint32_t* untigLen, 
long long minLongUntig, long long maxShortUntig, float ShortUntigRate, long long mainLen)
{
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    uint32_t Len, endNode;
    if(untigLen == NULL)
    {
        ///Len = get_long_tip_length(sg, u, begNode, &endNode, NULL);
        if(get_unitig(sg, ug, begNode, &endNode, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, NULL) == LOOP)
        {
            ///the length of LOOP is infinite
            return 1;
        }
        Len = nodeLen;
    }
    else
    {
        Len = (*untigLen);
    }
    
    
    if(Len == 0) return 0;
    if(Len < (ShortUntigRate*mainLen)) return 0;
    if(Len >= maxShortUntig) return 1;
    if(Len >= minLongUntig && Len >= (ShortUntigRate*mainLen)) return 1;
    return 0;
}


long long check_if_diploid_untigs(asg_t *nsg, asg_t *read_sg, uint32_t v_0, uint32_t v_1, 
ma_utg_v* ut_v,  ma_hit_t_alloc* reverse_sources, long long min_edge_length, 
float hap_rate, buf_t* b_0, buf_t* b_1, int if_drop, R_to_U* ruIndex)
{
    uint32_t vEnd;
    b_0->b.n = b_1->b.n = 0;
    if(get_long_tip_length(nsg, ut_v, v_0, &vEnd, b_0) == 0) return -1;
    if(get_long_tip_length(nsg, ut_v, v_1, &vEnd, b_1) == 0) return -1;
    uint32_t l_0 = 0, l_1 = 0, i;
    for (i = 0; i < b_0->b.n; i++)
    {
        l_0 += ut_v->a[(b_0->b.a[i])>>1].n;
    }

    for (i = 0; i < b_1->b.n; i++)
    {
        l_1 += ut_v->a[(b_1->b.a[i])>>1].n;
    }
    
    if((long long)(l_0) <= min_edge_length || (long long)(l_1) <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_max;
    buf_t* b_min;
    
    if(l_0 <= l_1)
    {
        b_min = b_0;
        b_max = b_1;
    }
    else
    {
        b_min = b_1;
        b_max = b_0;
    }


    

    uint32_t max_count = 0;
    uint32_t min_count = 0;
    ma_utg_t* node_a;
    uint32_t i_a_1, i_a_2, untigID_a, qn, tn, j;
    ma_utg_t* node_b;
    uint32_t i_b_1, i_b_2, untigID_b;
    uint32_t is_Unitig;
    for (i_a_1 = 0; i_a_1 < b_min->b.n; i_a_1++)
    {
        untigID_a = b_min->b.a[i_a_1]>>1;
        node_a = &(ut_v->a[untigID_a]);
        for (i_a_2 = 0; i_a_2 < node_a->n; i_a_2++)  // iterate over the unitig's reads
        {
            qn = node_a->a[i_a_2]>>33;

            for (j = 0; j < reverse_sources[qn].length; j++)
            {
                tn = Get_tn(reverse_sources[qn].buffer[j]);
                /****************************may have bugs********************************/
                ///here is bug
                // if(read_sg->seq[tn].del == 1 || (if_drop == 1 && read_sg->seq[tn].c == ALTER_LABLE)) 
                // {
                //     continue;
                // }
                if(read_sg->seq[tn].del == 1)
                {
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
                }
                /****************************may have bugs********************************/
                min_count++;



                for (i_b_1 = 0; i_b_1 < b_max->b.n; i_b_1++)
                {
                    untigID_b = b_max->b.a[i_b_1]>>1;
                    node_b = &(ut_v->a[untigID_b]);
                    for (i_b_2 = 0; i_b_2 < node_b->n; i_b_2++)
                    {
                        if(tn == (node_b->a[i_b_2]>>33))
                        {
                            max_count++;
                            goto end_reverse;
                        }
                    }
                }
                end_reverse:;
            }
        }
    }
    


    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count > min_count*hap_rate) return 1;
    return 0;
}

long long check_if_diploid_untigs_complex(asg_t *nsg, asg_t *read_sg, uint32_t v_0, uint32_t v_1, 
ma_utg_v* ut_v,  ma_hit_t_alloc* reverse_sources, long long min_edge_length, long long stops_threshold,
float hap_rate, buf_t* b_0, buf_t* b_1, int if_drop, R_to_U* ruIndex)
{
    uint32_t vEnd;
    b_0->b.n = b_1->b.n = 0;

    
    if(get_long_tip_length_stops(nsg, ut_v, v_0, &vEnd, b_0, stops_threshold) == 0) return -1;
    if(get_long_tip_length_stops(nsg, ut_v, v_1, &vEnd, b_1, stops_threshold) == 0) return -1;


    uint32_t l_0 = 0, l_1 = 0;
    long long i, j;
    i = b_0->b.n; i--;
    j = b_1->b.n; j--;
    while (i>=0 && j >=0)
    {
        if(b_0->b.a[i] == b_1->b.a[j])
        {
            i--;
            j--;
        }
        else
        {
            break;
        }
    }
    b_0->b.n = i+1;
    b_1->b.n = j+1;


    l_0 = 0; l_1 = 0;
    for (i = 0; i < (long long)b_0->b.n; i++)
    {
        l_0 += ut_v->a[(b_0->b.a[i])>>1].n;
    }

    for (i = 0; i < (long long)b_1->b.n; i++)
    {
        l_1 += ut_v->a[(b_1->b.a[i])>>1].n;
    }
    
    if((long long)(l_0) <= min_edge_length || (long long)(l_1) <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_max;
    buf_t* b_min;
    
    if(l_0 <= l_1)
    {
        b_min = b_0;
        b_max = b_1;
    }
    else
    {
        b_min = b_1;
        b_max = b_0;
    }


    

    uint32_t max_count = 0;
    uint32_t min_count = 0;
    ma_utg_t* node_a;
    uint32_t i_a_1, i_a_2, untigID_a, qn, tn;
    ma_utg_t* node_b;
    uint32_t i_b_1, i_b_2, untigID_b;
    uint32_t is_Unitig;
    for (i_a_1 = 0; i_a_1 < b_min->b.n; i_a_1++)
    {
        untigID_a = b_min->b.a[i_a_1]>>1;
        node_a = &(ut_v->a[untigID_a]);
        for (i_a_2 = 0; i_a_2 < node_a->n; i_a_2++)
        {
            qn = node_a->a[i_a_2]>>33;

            for (j = 0; j < (long long)reverse_sources[qn].length; j++)
            {
                tn = Get_tn(reverse_sources[qn].buffer[j]);
                ///here is bug
                /****************************may have bugs********************************/
                // if(read_sg->seq[tn].del == 1 || (if_drop == 1 && read_sg->seq[tn].c == ALTER_LABLE)) 
                // {
                //     continue;
                // }
                if(read_sg->seq[tn].del == 1)
                {
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_sg->seq[tn].del == 1) continue;
                }
                /****************************may have bugs********************************/
                min_count++;



                for (i_b_1 = 0; i_b_1 < b_max->b.n; i_b_1++)
                {
                    untigID_b = b_max->b.a[i_b_1]>>1;
                    node_b = &(ut_v->a[untigID_b]);
                    for (i_b_2 = 0; i_b_2 < node_b->n; i_b_2++)
                    {
                        if(tn == (node_b->a[i_b_2]>>33))
                        {
                            max_count++;
                            goto end_reverse;
                        }
                    }
                }
                end_reverse:;
            }
        }
    }
    


    if(min_count == 0) return -1;
    if(max_count == 0) return 0;
    if(max_count > min_count*hap_rate) return 1;
    return 0;
}



long long check_if_diploid_aggressive(uint32_t v1, uint32_t v2, asg_t *g, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    uint32_t convex1, convex2;
    long long l1, l2;

    b_0.b.n = 0;
    b_1.b.n = 0;
    ///uint32_t flag1 = detect_single_path(g, v1, &convex1, &l1, &b_0);
    uint32_t flag1 = detect_single_path_with_dels(g, v1, &convex1, &l1, &b_0);
    
    ///uint32_t flag2 = detect_single_path(g, v2, &convex2, &l2, &b_1);
    uint32_t flag2 = detect_single_path_with_dels(g, v2, &convex2, &l2, &b_1);
    
    if(flag1 == LOOP || flag2 == LOOP)
    {
        return -1;
    }

    if(flag1 != END_TIPS && flag1 != LONG_TIPS)
    {
        l1--;
        b_0.b.n--;
    }

    if(flag2 != END_TIPS && flag2 != LONG_TIPS)
    {
        l2--;
        b_1.b.n--;
    }


    if(l1 <= min_edge_length || l2 <= min_edge_length)
    {
        return -1;
    }

    buf_t* b_min;
    buf_t* b_max;
    if(l1<=l2)
    {
        b_min = &b_0;
        b_max = &b_1;
    }
    else
    {
        b_min = &b_1;
        b_max = &b_0;
    }

    long long i, j, k;
    double max_count = 0;
    double min_count = 0;
    uint32_t qn, tn;
    for (i = 0; i < (long long)b_min->b.n; i++)
    {
        qn = b_min->b.a[i];
        for (j = 0; j < (long long)reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(g->seq[tn].del == 1) continue;
            min_count++;
            for (k = 0; k < (long long)b_max->b.n; k++)
            {
                if(b_max->b.a[k]==tn)
                {
                    max_count++;
                    break;
                }
            }
        }
    }



    free(b_0.b.a);
    free(b_1.b.a);

    if(min_count == 0) return -1;
    if(max_count == 0) return 0;

    return 1;
    /**
    if(max_count/min_count>0.3) return 1;
    return 0;
    **/

}

int asg_arc_del_too_short_overlaps(asg_t *g, long long dropLen, float drop_ratio, 
ma_hit_t_alloc* reverse_sources, long long min_edge_length, R_to_U* ruIndex)
{
    double startTime = Get_T();

	uint32_t v, v_max, v_maxLen, n_vtx = g->n_seq * 2, n_short = 0;
    long long drop_ratio_Len = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;
        if (g->seq_vis[v] != 0) continue;

        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv < 2 || g->seq[v>>1].del) continue;

        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;
        v_max = (uint32_t)-1;

        for (i = 0, n_arc = 0; i < nv; i++)  // iterates from longest ovlp to the shortest (since sorting v_lv is sorting lu)
        {
            if (!av[i].del)
            {
                if(v_max == (uint32_t)-1)  // 1st
                {
                    v_max = av[i].v;
                    v_maxLen = av[i].ol;
                    if(v_maxLen < dropLen) break;
                    drop_ratio_Len = v_maxLen * drop_ratio;
                    if(dropLen < drop_ratio_Len)
                    {
                        drop_ratio_Len = dropLen;
                    }
                }
                else if(av[i].ol < drop_ratio_Len && 
                check_if_diploid(v_max, av[i].v, g, reverse_sources, min_edge_length, ruIndex) != 1)
                {
                    av[i].ol = 1;
                    asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
                    ++n_short;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
	
    fprintf(stderr, "[M::%s] removed %d short overlaps, used %.1fs\n", 
                    __func__, n_short, Get_T()-startTime);

	return n_short;
}



int asg_arc_del_short_diploid_unclean_exact(asg_t *g, float drop_ratio, ma_hit_t_alloc* sources)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
	for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del) continue;

		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, nv = asg_arc_n(g, v);
		///if there is just one overlap, do nothing
		if (nv < 2) continue;
		///keep the longest one
        for (i = 1; i < nv; i++)
        {
            ///if it is an inexact overlap
            if(av[i].el == 0 && 
              sources[v>>1].is_fully_corrected == 1&& 
              sources[(av[i].v>>1)].is_fully_corrected == 1)
            {
                av[i].del = 1;
                ++n_short;
            }
        }

	}

	if (n_short) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d inexact overlaps\n", __func__, n_short);
	return n_short;
}



///check if v has only one branch
static uint32_t asg_check_unambi1(asg_t *g, uint32_t v)
{
	asg_arc_t *av = asg_arc_a(g, v);
	uint32_t i, nv = asg_arc_n(g, v);
	uint32_t k = nv, kv;
	for (i = 0, kv = 0; i < nv; ++i)
		if (!av[i].del) ++kv, k = i;
	if (kv != 1) return (uint32_t)-1;
	return av[k].v;
}
///to see if it is a long tip
static int asg_topocut_aux(asg_t *g, uint32_t v, int max_ext)
{
	int32_t n_ext;
	for (n_ext = 1; n_ext < max_ext && v != (uint32_t)-1; ++n_ext) {
		if (asg_check_unambi1(g, v^1) == (uint32_t)-1) {
			--n_ext;
			break;
		}
		v = asg_check_unambi1(g, v);
	}
    
	return n_ext;
}

// delete short arcs
///for best graph?
int asg_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio, int max_ext, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold,
int if_skip_bubble, int if_drop, int if_check_hap, R_to_U* ruIndex)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(if_skip_bubble && g->seq_vis[v] != 0) continue;
        if(if_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0, ow_max_i = 0;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) ov_max = av[i].ol, ov_max_i = i; 
			++kv;
		}
		if (kv >= 2 && a->ol > ov_max * drop_ratio) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) ow_max = aw[i].ol, ow_max_i = i;
			++kw;
		}
		if (kw >= 2 && a->ol > ow_max * drop_ratio) continue;
		///if (kv == 1 && kw == 1) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        ///kv and kw is the avialiable 
		if (kv > 1 && kw > 1) {
			if (a->ol < ov_max * drop_ratio && a->ol < ow_max * drop_ratio)
				to_del = 1;
            if(if_check_hap == 1 && to_del == 1)
            {
                
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1 
                   || 
                   check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
        
		} else if (kw == 1) {
			if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            ///kv > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
		} else if (kv == 1) {
			if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            ///kw > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
		}
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %lld short overlaps, used %.1f\n", 
                    __func__, n_cut, Get_T()-startTime);
	return n_cut;
}

int unitig_arc_del_short_diploid_by_length_topo(asg_t *g, ma_ug_t *ug, float drop_ratio, 
int max_ext, ma_hit_t_alloc* reverse_sources, int if_skip_bubble, int if_drop)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2, convex;
    long long n_cut = 0, tmp, nodeLen, max_stop_nodeLen, max_stop_baseLen;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(if_skip_bubble && g->seq_vis[v] != 0) continue;
        if(if_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        if(get_real_length(g, v, NULL) < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) ov_max = av[i].ol; 
			++kv;
		}
		if (kv >= 2 && a->ol > ov_max * drop_ratio) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) ow_max = aw[i].ol;
			++kw;
		}
		if (kw >= 2 && a->ol > ow_max * drop_ratio) continue;
		///if (kv == 1 && kw == 1) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        ///kv and kw is the avialiable 
		if (kv > 1 && kw > 1) {
			if (a->ol < ov_max * drop_ratio && a->ol < ow_max * drop_ratio)
            {
                to_del = 1;
            }
		} else if (kw == 1) {
            get_unitig(g, ug, w^1, &convex, &nodeLen, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 1, NULL);
            if(nodeLen < max_ext) to_del = 1;
            ///if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
		} else if (kv == 1) {
            get_unitig(g, ug, v^1, &convex, &nodeLen, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 1, NULL);
            if(nodeLen < max_ext) to_del = 1;
			///if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
		}
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld short overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}


int unitig_arc_del_short_diploid_by_length(asg_t *g, float drop_ratio)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq[v>>1].c == ALTER_LABLE || g->seq[v>>1].del) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, v = (a->ul)>>32;
		uint32_t nv = asg_arc_n(g, v), kv;
		uint32_t ov_max = 0;
		asg_arc_t *av = NULL;
		///nv must be >= 2
		if (nv <= 1) continue;
        av = asg_arc_a(g, v);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) ov_max = av[i].ol; 
			++kv;
		}
        if (kv <= 1) continue;
		if (kv >= 2 && a->ol > ov_max * drop_ratio) continue;
        a->del = 1;
        asg_arc_del(g, a->v^1, av->ul>>32^1, 1);
        ++n_cut; 
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld short overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}


uint8_t get_tip_trio_infor(asg_t *sg, uint32_t begNode)
{
    uint32_t v = begNode, w;
    uint32_t kv;
    uint32_t eLen = 0, uLen = 0;
    uint32_t father_occ = 0, mother_occ = 0, ambigious_occ = 0;

    while (1)
    {
        kv = get_real_length(sg, v, NULL);
        eLen++;

        if(R_INF.trio_flag[v>>1]==FATHER)
        {
            father_occ++;
        }
        else if(R_INF.trio_flag[v>>1]==MOTHER)
        {
            mother_occ++;
        }
        else if((R_INF.trio_flag[v>>1]==AMBIGU) || (R_INF.trio_flag[v>>1]==DROP))
        {
            ambigious_occ++;
        }
        
        if(kv!=1) break;
        ///kv must be 1 here
        kv = get_real_length(sg, v, &w);
        if(get_real_length(sg, w^1, NULL)!=1) break;
        v = w;
        if(v == begNode) break;
    }

    uLen = eLen;
    eLen = father_occ + mother_occ;
    if(eLen == 0) return AMBIGU;
    if(father_occ >= mother_occ)
    {
        if((father_occ > TRIO_THRES*eLen) && (father_occ >= DOUBLE_CHECK_THRES*uLen)) return FATHER;
    }
    else
    {
        if((mother_occ > TRIO_THRES*eLen) && (mother_occ >= DOUBLE_CHECK_THRES*uLen)) return MOTHER;
    }
    return AMBIGU;
}

int asg_arc_del_short_diploid_by_length_trio(asg_t *g, float drop_ratio, int max_ext, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold,
int if_skip_bubble, int if_drop, int if_check_hap, R_to_U* ruIndex)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(if_skip_bubble && g->seq_vis[v] != 0) continue;
        if(if_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        uint64_t i;
        for (i = 0; i < nv; ++i)
        {
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0, ow_max_i = 0, trio_flag, non_trio_flag;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);
        kv = get_real_length(g, v, NULL);
        kw = get_real_length(g, w, NULL);
        if (kv <= 1 && kw <= 1) continue;
        trio_flag = get_tip_trio_infor(g, v^1);
        non_trio_flag = (uint32_t)-1;
        if(trio_flag == FATHER) non_trio_flag = MOTHER;
        if(trio_flag == MOTHER) non_trio_flag = FATHER;
        
        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
            kv++; 
            if(get_tip_trio_infor(g, av[i].v) == non_trio_flag) continue;
			if (ov_max < av[i].ol) ov_max = av[i].ol, ov_max_i = i;
            ///kv++; 
		}
		if (kv >= 2 && a->ol > ov_max * drop_ratio) continue;

		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
            kw++;
            if (get_tip_trio_infor(g, aw[i].v) == non_trio_flag) continue;
			if (ow_max < aw[i].ol) ow_max = aw[i].ol, ow_max_i = i;
            ///kw++;
		}
		if (kw >= 2 && a->ol > ow_max * drop_ratio) continue;

        if (kv <= 1 && kw <= 1) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        ///kv and kw is the avialiable 
		if (kv > 1 && kw > 1) {
			if (a->ol < ov_max * drop_ratio && a->ol < ow_max * drop_ratio)
				to_del = 1;
            if(if_check_hap == 1 && to_del == 1)
            {
                
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1 
                   || 
                   check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
        
		} else if (kw == 1) {
			if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            ///kv > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(av[ov_max_i].v, w^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
		} else if (kv == 1) {
			if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            ///kw > 1
            if(if_check_hap == 1 && to_del == 1)
            {
                if(check_if_diploid_primary_complex(aw[ow_max_i].v, v^1, g, reverse_sources, 
                miniedgeLen, stops_threshold, if_drop, ruIndex) == 1)
                {
                    to_del = 0;
                }
            }
		}
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld short overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}

int asg_arc_del_short_false_link(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));

	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del_l, to_del_r;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del_l = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                ///kv_push(uint32_t, b_r, aw[t].v);
            }

            if(av[i].ol < min_edge * drop_ratio) to_del_l++;
        }












        /****************************may have bugs********************************/
        if(to_del_l != kv)
        {
            b_f.n = 0;
            b_r.n = 0;
            to_del_l = 0;

            for (i = 0; i < nv; i++)
            {
                if (av[i].del) continue;

                w = av[i].v^1;
                nw = asg_arc_n(g, w);
                if(nw < 2) break;
                kw = get_real_length(g, w, NULL);
                if(kw < 2) break;

                kv_push(uint32_t, b_f, av[i].v);
                kv_push(uint32_t, b_r, w);

                aw = asg_arc_a(g, w);
                min_edge = (u_int32_t)-1;
                for (t = 0; t < nw; t++)
                {
                    if(aw[t].del) continue;
                    if((aw[t].v>>1) == (v>>1)) continue;
                    if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                }

                if(av[i].ol < min_edge * o_drop_ratio) to_del_l++;
            }

            if(to_del_l == kv)
            {
                ///forward
                to_del_l = 1;
                for (i = 1; i < b_f.n; i++)
                {
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


                    w = b_r.a[0];
                    kw = get_real_length(g, w, NULL);
                    if(kw != 2) goto terminal;
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (t = 0; t < nw; t++)
                    {
                        if(aw[t].del) continue;
                        if((aw[t].v>>1) == (v>>1)) continue;
                        w0 = aw[t].v;
                    }
                    to_del_l = 1;


                    for (i = 1; i < b_r.n; i++)
                    {
                        w = b_r.a[i];
                        kw = get_real_length(g, w, NULL);
                        if(kw != 2) goto terminal;
                        aw = asg_arc_a(g, w);
                        nw = asg_arc_n(g, w);
                        for (t = 0; t < nw; t++)
                        {
                            if(aw[t].del) continue;
                            if((aw[t].v>>1) == (v>>1)) continue;
                            w1 = aw[t].v;
                        }

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
                        {
                            to_del_l++;
                        }
                    }

                }
            }
        }

        terminal:
        /****************************may have bugs********************************/


        if(to_del_l != kv) continue;





        uint32_t convex1;
        long long l1;




        ////forward bubble
        to_del_l = 0;
        for (i = 0; i < b_f.n; i++)
        {
            if(b_f.a[i] == b_f.a[0])
            {
                to_del_l = 1;
            } 
            else
            {
                to_del_l = 0;
                break;
            }
        }
        //check the length
        if(to_del_l == 0 && asg_bub_end_finder_with_del_advance(g, 
        b_f.a, b_f.n, max_dist, &bub, 0, (u_int32_t)-1, &sink)==1)
        {
            to_del_l = 1;
        }
        if(to_del_l == 0 && detect_mul_bubble_end_with_bubbles(g, b_f.a, b_f.n, &convex1, &l1, NULL))
        {
            to_del_l = 1;
        }

        ////backward bubble
        to_del_r = 0;
        for (i = 0; i < b_r.n; i++)
        {
            if(b_r.a[i] == b_r.a[0])
            {
                to_del_r = 1;
            } 
            else
            {
                to_del_r = 0;
                break;
            }
        }
        if(to_del_r == 0 && asg_bub_end_finder_with_del_advance
        (g, b_r.a, b_r.n, max_dist, &bub, 1, v^1, &sink)==1)
        {
            to_del_r = 1;
        }
        if(to_del_r == 0 && detect_mul_bubble_end_with_bubbles(g, b_r.a, b_r.n, &convex1, &l1, NULL))
        {
            to_del_r = 1;
        }

        
    



		if (to_del_l && to_del_r)
        {
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %u false overlaps, used %.1fs\n", 
                    __func__, n_cut, Get_T()-startTime);
	return n_cut;
}


int asg_arc_del_short_false_link_primary(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq[v>>1].c == ALTER_LABLE) continue;

        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del_l, to_del_r;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del_l = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                ///kv_push(uint32_t, b_r, aw[t].v);
            }

            if(av[i].ol < min_edge * drop_ratio) to_del_l++;
        }












        /****************************may have bugs********************************/
        if(to_del_l != kv)
        {
            b_f.n = 0;
            b_r.n = 0;
            to_del_l = 0;

            for (i = 0; i < nv; i++)
            {
                if (av[i].del) continue;

                w = av[i].v^1;
                nw = asg_arc_n(g, w);
                if(nw < 2) break;
                kw = get_real_length(g, w, NULL);
                if(kw < 2) break;

                kv_push(uint32_t, b_f, av[i].v);
                kv_push(uint32_t, b_r, w);

                aw = asg_arc_a(g, w);
                min_edge = (u_int32_t)-1;
                for (t = 0; t < nw; t++)
                {
                    if(aw[t].del) continue;
                    if((aw[t].v>>1) == (v>>1)) continue;
                    if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                }

                if(av[i].ol < min_edge * o_drop_ratio) to_del_l++;
            }

            if(to_del_l == kv)
            {
                ///forward
                to_del_l = 1;
                for (i = 1; i < b_f.n; i++)
                {
                    /****************************may have bugs********************************/
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {/****************************may have bugs********************************/
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


                    w = b_r.a[0];
                    kw = get_real_length(g, w, NULL);
                    if(kw != 2) goto terminal;
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (t = 0; t < nw; t++)
                    {
                        if(aw[t].del) continue;
                        if((aw[t].v>>1) == (v>>1)) continue;
                        w0 = aw[t].v;
                    }
                    to_del_l = 1;


                    for (i = 1; i < b_r.n; i++)
                    {
                        w = b_r.a[i];
                        kw = get_real_length(g, w, NULL);
                        if(kw != 2) goto terminal;
                        aw = asg_arc_a(g, w);
                        nw = asg_arc_n(g, w);
                        for (t = 0; t < nw; t++)
                        {
                            if(aw[t].del) continue;
                            if((aw[t].v>>1) == (v>>1)) continue;
                            w1 = aw[t].v;
                        }
                        /****************************may have bugs********************************/
                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
                        {/****************************may have bugs********************************/
                            to_del_l++;
                        }
                    }

                }
            }
        }

        terminal:
        /****************************may have bugs********************************/


        if(to_del_l != kv) continue;





        uint32_t convex1;
        long long l1;




        ////forward bubble
        to_del_l = 0;
        for (i = 0; i < b_f.n; i++)
        {
            if(b_f.a[i] == b_f.a[0])
            {
                to_del_l = 1;
            } 
            else
            {
                to_del_l = 0;
                break;
            }
        }
        //check the length
        if(to_del_l == 0 && asg_bub_end_finder_with_del_advance(g, 
        b_f.a, b_f.n, max_dist, &bub, 0, (u_int32_t)-1, &sink)==1)
        {
            to_del_l = 1;
        }
        if(to_del_l == 0 && detect_mul_bubble_end_with_bubbles(g, b_f.a, b_f.n, &convex1, &l1, NULL))
        {
            to_del_l = 1;
        }

        ////backward bubble
        to_del_r = 0;
        for (i = 0; i < b_r.n; i++)
        {
            if(b_r.a[i] == b_r.a[0])
            {
                to_del_r = 1;
            } 
            else
            {
                to_del_r = 0;
                break;
            }
        }
        if(to_del_r == 0 && asg_bub_end_finder_with_del_advance
        (g, b_r.a, b_r.n, max_dist, &bub, 1, v^1, &sink)==1)
        {
            to_del_r = 1;
        }
        if(to_del_r == 0 && detect_mul_bubble_end_with_bubbles(g, b_r.a, b_r.n, &convex1, &l1, NULL))
        {
            to_del_r = 1;
        }

        
    



		if (to_del_l && to_del_r)
        {
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %u false overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}

int asg_arc_del_short_false_link_advance(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;
    uint32_t sink;

    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    uint64_t k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del_l, to_del_r;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del_l = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                ///kv_push(uint32_t, b_r, aw[t].v);
            }

            if(av[i].ol < min_edge * drop_ratio) to_del_l++;
        }












        /****************************may have bugs********************************/
        if(to_del_l != kv)
        {
            b_f.n = 0;
            b_r.n = 0;
            to_del_l = 0;

            for (i = 0; i < nv; i++)
            {
                if (av[i].del) continue;

                w = av[i].v^1;
                nw = asg_arc_n(g, w);
                if(nw < 2) break;
                kw = get_real_length(g, w, NULL);
                if(kw < 2) break;

                kv_push(uint32_t, b_f, av[i].v);
                kv_push(uint32_t, b_r, w);

                aw = asg_arc_a(g, w);
                min_edge = (u_int32_t)-1;
                for (t = 0; t < nw; t++)
                {
                    if(aw[t].del) continue;
                    if((aw[t].v>>1) == (v>>1)) continue;
                    if(aw[t].ol < min_edge) min_edge = aw[t].ol;
                }

                if(av[i].ol < min_edge * o_drop_ratio) to_del_l++;
            }

            if(to_del_l == kv)
            {
                ///forward
                to_del_l = 1;
                for (i = 1; i < b_f.n; i++)
                {
                    if(check_if_diploid(b_f.a[0], b_f.a[i], g, reverse_sources, miniedgeLen, ruIndex) == 1)
                    {
                        to_del_l++;
                    }
                }

                ///backward
                if(to_del_l != kv && b_r.n >= 2)
                {
                    to_del_l = 0;
                    uint32_t w0 = 0, w1 = 0;


                    w = b_r.a[0];
                    kw = get_real_length(g, w, NULL);
                    if(kw != 2) goto terminal;
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (t = 0; t < nw; t++)
                    {
                        if(aw[t].del) continue;
                        if((aw[t].v>>1) == (v>>1)) continue;
                        w0 = aw[t].v;
                    }
                    to_del_l = 1;


                    for (i = 1; i < b_r.n; i++)
                    {
                        w = b_r.a[i];
                        kw = get_real_length(g, w, NULL);
                        if(kw != 2) goto terminal;
                        aw = asg_arc_a(g, w);
                        nw = asg_arc_n(g, w);
                        for (t = 0; t < nw; t++)
                        {
                            if(aw[t].del) continue;
                            if((aw[t].v>>1) == (v>>1)) continue;
                            w1 = aw[t].v;
                        }

                        if(check_if_diploid(w0, w1, g, reverse_sources, miniedgeLen, ruIndex) == 1)
                        {
                            to_del_l++;
                        }
                    }

                }
            }
        }

        terminal:
        /****************************may have bugs********************************/


        if(to_del_l != kv) continue;





        uint32_t convex1;
        long long l1;




        ////forward bubble
        to_del_l = 0;
        for (i = 0; i < b_f.n; i++)
        {
            if(b_f.a[i] == b_f.a[0])
            {
                to_del_l = 1;
            } 
            else
            {
                to_del_l = 0;
                break;
            }
        }
        //check the length
        if(to_del_l == 0 && asg_bub_end_finder_with_del_advance(g, 
        b_f.a, b_f.n, max_dist, &bub, 0, (u_int32_t)-1, &sink)==1)
        {
            to_del_l = 1;
        }
        if(to_del_l == 0 && detect_mul_bubble_end_with_bubbles(g, b_f.a, b_f.n, &convex1, &l1, NULL))
        {
            to_del_l = 1;
        }

        ////backward bubble
        to_del_r = 0;
        for (i = 0; i < b_r.n; i++)
        {
            if(b_r.a[i] == b_r.a[0])
            {
                to_del_r = 1;
            } 
            else
            {
                to_del_r = 0;
                break;
            }
        }
        if(to_del_r == 0 && asg_bub_end_finder_with_del_advance
        (g, b_r.a, b_r.n, max_dist, &bub, 1, v^1, &sink)==1)
        {
            to_del_r = 1;
        }
        if(to_del_r == 0 && detect_mul_bubble_end_with_bubbles(g, b_r.a, b_r.n, &convex1, &l1, NULL))
        {
            to_del_r = 1;
        }

  


		if (to_del_l && to_del_r)
        {
            ///fprintf(stderr, "%.*s\n", Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
            for (i = 0; i < nv; ++i) 
            {
			    if (av[i].del) continue;
                
                ++n_cut;
                av[i].del = 1;
                asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
		    }

        }
    }
    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d false overlaps\n", __func__, n_cut);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
	return n_cut;
}


int asg_arc_del_complex_false_link(asg_t *g, float drop_ratio, float o_drop_ratio, int max_dist, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen)
{
    double startTime = Get_T();
    
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));


    kvec_t(uint32_t) b_f;
    memset(&b_f, 0, sizeof(b_f));

    kvec_t(uint32_t) b_r;
    memset(&b_r, 0, sizeof(b_r));


	uint32_t v, w, n_vtx = g->n_seq * 2, n_cut = 0;


    buf_t bub;
    if (!g->is_symm) asg_symm(g);
    memset(&bub, 0, sizeof(buf_t));
    bub.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);

            if(nv == 1 && asg_arc_n(g, v^1) == 1) continue;

            uint64_t t_ol = 0;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                t_ol += av[i].ol;
            }
            kv_push(uint64_t, b, (uint64_t)(t_ol << 32 | v));  
        }
	}


    radix_sort_arch64(b.a, b.a + b.n);

    uint32_t min_edge;
    
    
    uint64_t k, t;
    for (k = 0; k < b.n; k++)
    {
        ///v is the node
        v = (uint32_t)b.a[k];
        if (g->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(g, v), nw, to_del;
        if (nv < 2) continue;
        uint32_t kv = get_real_length(g, v, NULL), kw;
        if (kv < 2) continue;///V must have two out-nodes
        uint32_t i;
        asg_arc_t *av = asg_arc_a(g, v), *aw;

        b_f.n = 0;
        b_r.n = 0;
        to_del = 0;
        for (i = 0; i < nv; i++)
        {
            if (av[i].del) continue;
            
            w = av[i].v^1;
            nw = asg_arc_n(g, w);
            if(nw < 2) break;
            kw = get_real_length(g, w, NULL);
            if(kw < 2) break;

            kv_push(uint32_t, b_f, av[i].v);
            kv_push(uint32_t, b_r, w);

            aw = asg_arc_a(g, w);
            min_edge = (u_int32_t)-1;
            for (t = 0; t < nw; t++)
            {
                if(aw[t].del) continue;
                if((aw[t].v>>1) == (v>>1)) continue;
                if(aw[t].ol < min_edge) min_edge = aw[t].ol;
            }

            if(av[i].ol < min_edge * drop_ratio) to_del++;
        }

        if(to_del != kv) continue;

        for (i = 0; i < nv; ++i) 
        {
            if (av[i].del) continue;
            ++n_cut;
            av[i].del = 1;
            asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
        }
    }


    if(n_cut > 0)
    {
        for (v = 0; v < n_vtx; ++v) 
        {
            uint32_t nv = asg_arc_n(g, v);
            if (g->seq[v>>1].del)
            {
                continue;
            } 

            if(nv < 2)
            {
                continue;
            }


            if(asg_bub_finder_without_del_advance(g, v, max_dist, &bub) == 1)
            {
                uint32_t i;
                g->seq_vis[v] = 3;
                g->seq_vis[v^1] = 3;
                for (i = 0; i < bub.b.n; i++)
                {
                    g->seq_vis[bub.b.a[i]] = 3;
                    g->seq_vis[bub.b.a[i]^1] = 3;
                }
            }
        }


        for (v = 0; v < n_vtx; ++v) 
        {
            if(g->seq_vis[v] == 3) continue;

            uint32_t nv = asg_arc_n(g, v);
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t i;
            for (i = 0; i < nv; ++i) 
            {
                if (av[i].del && g->seq_vis[av[i].v] != 3)
                {
                    av[i].del = 0;
                    asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 0);
                } 
            }
        }
    }
    


    
    free(b.a); free(b_f.a); free(b_r.a);
    free(bub.a); free(bub.S.a); free(bub.T.a); free(bub.b.a); free(bub.e.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %d false overlaps, used %.1fs\n", __func__, n_cut, Get_T()-startTime);
	return n_cut;
}


int asg_arc_del_short_diploid_by_exact(asg_t *g, int max_ext, ma_hit_t_alloc* sources)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ol) 
            {
                ov_max = av[i].ol;
                ov_max_i = i;
            }
			++kv;
		}
		if (kv >= 2 && a->ol == ov_max) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ol) 
            {
                ow_max = aw[i].ol;
            }
			++kw;
		}
		if (kw >= 2 && a->ol == ow_max) continue;

		///if (kv == 1 && kw == 1) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;


        ///if this edge is an inexact edge
        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 1)
        {
            if (kv > 1 && kw > 1) {
				to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 0)
        {
            /****************************may have bugs********************************/
            ///if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected)
            /****************************may have bugs********************************/
            if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected == 1)
            {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
                } else if (kv == 1) {
                    if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
                }
            }
        }


        
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %lld inexact overlaps, used %.1fs\n", 
                    __func__, n_cut, Get_T()-startTime);
	return n_cut;
}

int asg_arc_del_short_diploid_by_exact_trio(asg_t *g, int max_ext, ma_hit_t_alloc* sources)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;
            long long i;
            for (i = 0; i < nv; ++i)
            {
                kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
            }
        }
	}

    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0, ov_max_i = 0, trio_flag, non_trio_flag;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);
        kv = get_real_length(g, v, NULL);
        kw = get_real_length(g, w, NULL);
        if (kv <= 1 && kw <= 1) continue;
        trio_flag = get_tip_trio_infor(g, v^1);
        non_trio_flag = (uint32_t)-1;
        if(trio_flag == FATHER) non_trio_flag = MOTHER;
        if(trio_flag == MOTHER) non_trio_flag = FATHER;


        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
            ++kv;
            if(get_tip_trio_infor(g, av[i].v) == non_trio_flag) continue;
			if (ov_max < av[i].ol) 
            {
                ov_max = av[i].ol;
                ov_max_i = i;
            }
			///++kv;
		}
		if (kv >= 2 && a->ol == ov_max) continue;


		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
            ++kw;
            if (get_tip_trio_infor(g, aw[i].v) == non_trio_flag) continue;
			if (ow_max < aw[i].ol) 
            {
                ow_max = aw[i].ol;
            }
			///++kw;
		}
		if (kw >= 2 && a->ol == ow_max) continue;
        if (kv <= 1 && kw <= 1) continue;


        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;


        ///if this edge is an inexact edge
        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 1)
        {
            if (kv > 1 && kw > 1) {
				to_del = 1;
            } else if (kw == 1) {
                if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
            } else if (kv == 1) {
                if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
            }
        }

        if(a->el == 0 && 
        sources[v>>1].is_fully_corrected == 1 && 
        sources[w>>1].is_fully_corrected == 0)
        {
            /****************************may have bugs********************************/
            ///if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected)
            /****************************may have bugs********************************/
            if(av[ov_max_i].el == 1 && sources[av[ov_max_i].v>>1].is_fully_corrected == 1)
            {
                if (kv > 1 && kw > 1) {
                    to_del = 1;
                } else if (kw == 1) {
                    if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
                } else if (kv == 1) {
                    if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
                }
            }
        }


        
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
        
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld inexact overlaps\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}


int asg_arc_del_short_diploi_by_suspect_edge(asg_t *g, int max_ext)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        ///if(g->seq_vis[v] == 0)
        {
            asg_arc_t *av = asg_arc_a(g, v);
            uint32_t nv = asg_arc_n(g, v);
            if (nv < 2) continue;

            long long i;
            for (i = 0; i < nv; ++i)
            {
                ///means there is a large indel at this edge
                if(av[i].no_l_indel == 0)
                {
                    kv_push(uint64_t, b, (uint64_t)((uint64_t)av[i].ol << 32 | (av - g->arc + i)));   
                }
            }
        }
	}


    radix_sort_arch64(b.a, b.a + b.n);


    uint64_t k;
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1, to_del = 0;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
		///nv must be >= 2
		if (nv == 1 && nw == 1) continue;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);

        ///calculate the longest edge for v and w
		for (i = 0, kv = 0; i < nv; ++i) 
        {
			if (av[i].del) continue;
			++kv;
		}

		for (i = 0, kw = 0; i < nw; ++i) 
        {
			if (aw[i].del) continue;
			++kw;
		}

		if (kv == 1 && kw == 1) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;

        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        
        if (kv > 1 && kw > 1) {
			to_del = 1;
        } else if (kw == 1) {
            if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del = 1;
        } else if (kv == 1) {
            if (asg_topocut_aux(g, v^1, max_ext) < max_ext) to_del = 1;
        }
        
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}
    
    fprintf(stderr, "[M::%s] removed %lld suspect overlaps, used %.1fs\n", 
                    __func__, n_cut, Get_T()-startTime);
	return n_cut;
}

int if_potential_false_node(ma_hit_t_alloc* paf, uint64_t rLen, ma_sub_t* max_left, ma_sub_t* max_right, int if_set)
{
    max_left->s = max_right->s = rLen;
    max_left->e = max_right->e = 0;

    long long j;
    uint32_t qs, qe;
    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;
        if(paf->buffer[j].el != 1) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);

        ///overlaps from left side
        if(qs == 0)
        {
            if(qs < max_left->s) max_left->s = qs;
            if(qe > max_left->e) max_left->e = qe;
        }

        ///overlaps from right side
        if(qe == rLen)
        {
            if(qs < max_right->s) max_right->s = qs;
            if(qe > max_right->e) max_right->e = qe;
        }

        ///note: if (qs == 0 && qe == rLen)
        ///this overlap would be added to both b_left and b_right
        ///that is what we want
    }

    if (max_left->e > max_right->s) return 0;

    long long new_left_e, new_right_s;
    new_left_e = max_left->e;
    new_right_s = max_right->s;

    for (j = 0; j < paf->length; j++)
    {
        if(paf->buffer[j].del) continue;
        if(paf->buffer[j].el != 1) continue;

        qs = Get_qs(paf->buffer[j]);
        qe = Get_qe(paf->buffer[j]);
        ///check contained overlaps
        if(qs != 0 && qe != rLen)
        {
            ///[qs, qe), [max_left.s, max_left.e)
            if(qs < max_left->e && qe > max_left->e)
            {
                ///if(qe > max_left->e) max_left->e = qe;
                if(qe > max_left->e && qe > new_left_e) new_left_e = qe;
            }

            ///[qs, qe), [max_right.s, max_right.e)
            if(qs < max_right->s && qe > max_right->s)
            {
                ///if(qs < max_right->s) max_right->s = qs;
                if(qs < max_right->s && qs < new_right_s) new_right_s = qs;
            }



            if(if_set)
            {
                max_left->e = new_left_e;
                max_right->s = new_right_s;
            }
        }
    }

    max_left->e = new_left_e;
    max_right->s = new_right_s;

    if (max_left->e > max_right->s) return 0;

    return 1;
}

//reomve edge between two chromesomes
//this node must be a single read
int asg_arc_del_false_node(asg_t *g, ma_hit_t_alloc* sources, int max_ext)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));
    ma_sub_t max_left, max_right;

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            if(asg_arc_n(g, v)!=1 || asg_arc_n(g, v^1)!=1)  // require node to be in a single path
            {
                continue;
            }

            if(asg_is_single_edge(g, asg_arc_a(g, v)[0].v, v>>1) < 2)
            {
                continue;
            }

            if(asg_is_single_edge(g, asg_arc_a(g, v^1)[0].v, (v^1)>>1) < 2)
            {
                continue;
            }

            if(asg_arc_a(g, v)[0].el == 1)
            {
                continue;
            }

            if(if_potential_false_node(&(sources[v>>1]), g->seq[v>>1].len, &max_left, &max_right, 1)==0)
            {
                continue;
            }

            asg_arc_t *av = asg_arc_a(g, v);
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[0].ol << 32 | (av - g->arc)));
        }
	}

    // each entry is: upper 32 bits overlap length, lower 32 bits which "group" does the target belong to
    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    ///here all edges are inexact matches 
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w  // ? this block is effectively calling get_real_length...
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			++kv;
		}

		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			++kw;
		}

		if (kv < 1 || kw < 2) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        

        uint32_t el_edges = 0;
        ///there should be at least two available edges in aw
        for (i = 0; i < nw; i++)
        {
            if (aw[i].del) continue;

            if(i != iw && aw[i].el == 1)
            {
                el_edges++;
            }
        }

        if(el_edges > 0 && av[iv].el == 0)
        {
            asg_seq_del(g, v>>1);
            ++n_cut;
        }
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %lld single nodes\n", __func__, n_cut);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
	return n_cut;
}
int asg_arc_del_false_node_meta(asg_t *g, int max_ext)
{
    double startTime = Get_T();
    kvec_t(uint64_t) b;
    memset(&b, 0, sizeof(b));
    ma_sub_t max_left, max_right;

	uint32_t v, n_vtx = g->n_seq * 2;
    long long n_cut = 0;
    
    for (v = 0; v < n_vtx; ++v) 
    {
        if(g->seq_vis[v] == 0)
        {
            if(asg_arc_n(g, v)!=1 || asg_arc_n(g, v^1)!=1)  // require node to be in a single path
            {
                continue;
            }

            if(asg_is_single_edge(g, asg_arc_a(g, v)[0].v, v>>1) < 2)
            {
                continue;
            }

            if(asg_is_single_edge(g, asg_arc_a(g, v^1)[0].v, (v^1)>>1) < 2)
            {
                continue;
            }

            if(asg_arc_a(g, v)[0].el == 1)
            {
                continue;
            }

            // if(if_potential_false_node(&(sources[v>>1]), g->seq[v>>1].len, &max_left, &max_right, 1)==0)
            // {
            //     continue;
            // }

            asg_arc_t *av = asg_arc_a(g, v);
            kv_push(uint64_t, b, (uint64_t)((uint64_t)av[0].ol << 32 | (av - g->arc)));
        }
	}

    // each entry is: upper 32 bits overlap length, lower 32 bits which "group" does the target belong to
    radix_sort_arch64(b.a, b.a + b.n);

    uint64_t k;
    ///here all edges are inexact matches 
    for (k = 0; k < b.n; k++)
    {
        
        asg_arc_t *a = &g->arc[(uint32_t)b.a[k]];
		///v is self id, w is the id of another end
		uint32_t i, iv, iw, v = (a->ul)>>32, w = a->v^1;
		uint32_t nv = asg_arc_n(g, v), nw = asg_arc_n(g, w), kv, kw;
		asg_arc_t *av, *aw;
        av = asg_arc_a(g, v);
		aw = asg_arc_a(g, w);


        ///calculate the longest edge for v and w  // ? this block is effectively calling get_real_length...
		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			++kv;
		}

		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			++kw;
		}

		if (kv < 1 || kw < 2) continue;

        ///to see which one is the current edge (from v and w)
		for (iv = 0; iv < nv; ++iv)
			if (av[iv].v == (w^1)) break;
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].v == (v^1)) break;
        ///if one edge has been deleted, it should be deleted in both direction
        if (av[iv].del && aw[iw].del) continue;

        

        uint32_t el_edges = 0;
        ///there should be at least two available edges in aw
        for (i = 0; i < nw; i++)
        {
            if (aw[i].del) continue;

            if(i != iw && aw[i].el == 1)
            {
                el_edges++;
            }
        }

        if(el_edges > 0 && av[iv].el == 0)
        {
            asg_seq_del(g, v>>1);
            ++n_cut;
        }
    }
    
    free(b.a);
	if (n_cut) 
    {
		asg_cleanup(g);
		asg_symm(g);
	}

    fprintf(stderr, "[M::%s] removed %lld single nodes, used %.1f\n", 
                    __func__, n_cut, Get_T()-startTime);
	return n_cut;
}

// #define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
// #define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

ma_ug_t *ma_ug_gen(asg_t *g)
{
    asg_cleanup(g);
	int32_t *mark;
	uint32_t i, v, n_vtx = g->n_seq * 2;
	kdq_t(uint64_t) *q;  // queue
	ma_ug_t *ug;

	ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
	ug->g = asg_init();
	mark = (int32_t*)calloc(n_vtx, 4);
    int flag_circ = 0;

	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
        flag_circ = 0;
		uint32_t w, x, l, start, end, len;
		ma_utg_t *p;
        ///what's the usage of mark array
        ///mark array is used to mark if this node has already been included in a contig
        /****************************may have bugs********************************/
		///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
        if (g->seq[v>>1].del || mark[v]) continue;
        if (arc_cnt(g, v) == 0 && arc_cnt(g, (v^1)) != 0) continue;
        /****************************may have bugs********************************/
		mark[v] = 1;
		q->count = 0, start = v, end = v^1, len = 0;
		// forward
		w = v;
		while (1) {
			/**
			 * w----->x
			 * w<-----x
			 * that means the only suffix of w is x, and the only prefix of x is w
			 **/
			if (arc_cnt(g, w) != 1) break;
			x = arc_first(g, w).v; // w->x
			if (arc_cnt(g, x^1) != 1) break;

			/**
			 * another direction of w would be marked as used (since w has been used)
			**/
			mark[x] = mark[w^1] = 1;
			///l is the edge length, instead of overlap length
            ///note: edge length is different with overlap length
			l = asg_arc_len(arc_first(g, w));
			kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
			end = x^1, len += l;
			if (x == v) {flag_circ = 1; break;}  // keep current w, which will be the end vertex (use w^1)
            w = x;
		}
		if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
			///length of seq, instead of edge
			l = g->seq[end>>1].len;
			kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
			len += l;
		} else { // circular unitig
			start = v; end = w^1; /*flag_circ = 1*/; // start = end = UINT32_MAX;
			goto add_unitig; // then it is not necessary to do the backward
		}
		// backward
		x = v;
		while (1) { // similar to forward but not the same
			if (arc_cnt(g, x^1) != 1) break;
			w = arc_first(g, x^1).v ^ 1; // w->x
			if (arc_cnt(g, w) != 1) break;
			mark[x] = mark[w^1] = 1;
			l = asg_arc_len(arc_first(g, w));
			///w is the seq id + direction, l is the length of edge
			///push element to the front of a queue
			kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);

            // fprintf(stderr, "uId: %u, >%.*s (%u)\n", 
            // ug->u.n, (int)Get_NAME_LENGTH((R_INF), w>>1), Get_NAME((R_INF), w>>1), w>>1);

			start = w, len += l;
			x = w;
		}
add_unitig:
		// if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        mark[start] = mark[end] = 1;
		kv_pushp(ma_utg_t, ug->u, &p);
		p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = /*(start == UINT32_MAX)*/ flag_circ;
		p->m = p->n;
		kv_roundup32(p->m);
		p->a = (uint64_t*)malloc(8 * p->m);
		//all elements are saved here
		for (i = 0; i < kdq_size(q); ++i)
			p->a[i] = kdq_at(q, i);
	}
	kdq_destroy(uint64_t, q);

	// add arcs between unitigs; reusing mark for a different purpose
	//ug saves all unitigs
	for (v = 0; v < n_vtx; ++v) mark[v] = -1;

    //mark all start nodes and end nodes of all unitigs
	for (i = 0; i < ug->u.n; ++i) {
		if (ug->u.a[i].circ) continue;
		mark[ug->u.a[i].start] = i<<1 | 0;
		mark[ug->u.a[i].end] = i<<1 | 1;
	}

	//scan all edges
	for (i = 0; i < g->n_arc; ++i) {
		asg_arc_t *p = &g->arc[i];
		if (p->del) continue;
		///to connect two unitigs, we need to connect the end of unitig x to the start of unitig y
		///so we need to ^1 to get the reverse direction of (x's end)?
        ///>=0 means this node is a start/end node of an unitig
        ///means this node is a intersaction node
		if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
			asg_arc_t *q;
			uint32_t u = mark[p->ul>>32^1]^1;
			int l = ug->u.a[u>>1].len - p->ol;
			if (l < 0) l = 1;
			q = asg_arc_pushp(ug->g);
			q->ol = p->ol, q->del = 0;
			q->ul = (uint64_t)u<<32 | l;
			q->v = mark[p->v];
		}
	}

    // add arcs for circular utg
    // (the same as above)
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;  // recycle `mark`
	for (i = 0; i < ug->u.n; ++i) {  // collect start/end
		if (!ug->u.a[i].circ) continue;
		mark[ug->u.a[i].start] = i<<1 | 0;
		mark[ug->u.a[i].end] = i<<1 | 1;
	}
    for (i = 0; i < g->n_arc; ++i) {  // find edges
		asg_arc_t *p = &g->arc[i];
		if (p->del) continue;
		if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
			asg_arc_t *q;
			uint32_t u = mark[p->ul>>32^1]^1;
			int l = ug->u.a[u>>1].len - p->ol;
			if (l < 0) l = 1;
			q = asg_arc_pushp(ug->g);
			q->ol = p->ol, q->del = 0;
			q->ul = (uint64_t)u<<32 | l;
			q->v = mark[p->v];
		}
	}

	for (i = 0; i < ug->u.n; ++i)
		asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
	asg_cleanup(ug->g);
	free(mark);
	return ug;
}







ma_ug_t *ma_ug_gen_primary(asg_t *g, uint8_t flag)
{   // in ha used by
    //     label_tangles
    //     recover_utg_by_coverage
    //     output_contig_graph_primary_pre (p_utg)
    //     output_contig_graph_primary (p_ctg)
    //     output_contig_graph_alternative
    asg_cleanup(g);
	int32_t *mark;
	uint32_t i, v, n_vtx = g->n_seq * 2;
	kdq_t(uint64_t) *q;  // queue
	ma_ug_t *ug;

	ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
	ug->g = asg_init();
	mark = (int32_t*)calloc(n_vtx, 4);

    int flag_circ = 0;

    ///for each untig, all node have the same direction
    ///and all node except the last one just have one edge
    ///the last one may have multiple edges
	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
        flag_circ = 0;
		uint32_t w, x, l, start, end, len;
		ma_utg_t *p;
        ///what's the usage of mark array
        ///mark array is used to mark if this node has already been included in a contig
        /****************************may have hap bugs********************************/
		///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v] || g->seq[v>>1].c != flag) continue;
        ///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v] || (g->seq[v>>1].c & flag) == 0) continue;
        ///if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
        if (g->seq[v>>1].del || mark[v]) continue;
        if (arc_cnt(g, v) == 0 && arc_cnt(g, (v^1)) != 0) continue;
        if(flag == PRIMARY_LABLE && g->seq[v>>1].c == ALTER_LABLE) continue;
        if(flag == ALTER_LABLE && g->seq[v>>1].c != ALTER_LABLE) continue;
        /****************************may have hap bugs********************************/
		mark[v] = 1;
		q->count = 0, start = v, end = v^1, len = 0;
		// forward
		w = v;
		while (1) {
			/**
			 * w----->x
			 * w<-----x
			 * that means the only suffix of w is x, and the only prefix of x is w
			 **/
			if (arc_cnt(g, w) != 1) break;
			x = arc_first(g, w).v; // w->x
			if (arc_cnt(g, x^1) != 1) break;

			/**
			 * another direction of w would be marked as used (since w has been used)
			**/
			mark[x] = mark[w^1] = 1;
			///l is the edge length, instead of overlap length
            ///note: edge length is different with overlap length
			l = asg_arc_len(arc_first(g, w));
			kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
			end = x^1, len += l;
			// w = x;
			// if (x == v) break;
            if (x == v) {flag_circ = 1; break;}  // keep current w, which will be the end vertex (use w^1)
            w = x;
		}
        ///kdq_size(q) == 0 means there is just one read
		if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
			///length of seq, instead of edge
			l = g->seq[end>>1].len;
			kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
			len += l;
		} else { // circular unitig
			start = v; end = w^1;  // start = end = UINT32_MAX;
			goto add_unitig; // then it is not necessary to do the backward
		}
		// backward
		x = v;
		while (1) { // similar to forward but not the same
			if (arc_cnt(g, x^1) != 1) break;
			w = arc_first(g, x^1).v ^ 1; // w->x
			if (arc_cnt(g, w) != 1) break;
			mark[x] = mark[w^1] = 1;
			l = asg_arc_len(arc_first(g, w));
			///w is the seq id + direction, l is the length of edge
			///push element to the front of a queue
			kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
			start = w, len += l;
			x = w;
		}
add_unitig:
		// if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        mark[start] = mark[end] = 1;
		kv_pushp(ma_utg_t, ug->u, &p);
		// p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
        p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = /*(start == UINT32_MAX)*/ flag_circ;
		p->m = p->n;
		kv_roundup32(p->m);
		p->a = (uint64_t*)malloc(8 * p->m);
		//all elements are saved here
		for (i = 0; i < kdq_size(q); ++i)
			p->a[i] = kdq_at(q, i);
	}
	kdq_destroy(uint64_t, q);

	// add arcs between unitigs; reusing mark for a different purpose
	//ug saves all unitigs
	for (v = 0; v < n_vtx; ++v) mark[v] = -1;

    //mark all start nodes and end nodes of all unitigs
    ///note ug->u.a[i].start == ug->u.a[i].a[0]
    ///but ug->u.a[i].end = ug->u.a[i].a[ug->u.a[i].n-1]^1
    ///ug->u.a[i].start has the same direction as other non-end node
    ///while ug->u.a[i].end is the only one with reverse direction
	for (i = 0; i < ug->u.n; ++i) {
		if (ug->u.a[i].circ) continue;
        ///i is the untig id
		mark[ug->u.a[i].start] = i<<1 | 0;
		mark[ug->u.a[i].end] = i<<1 | 1;
	}

	//scan all edges
	for (i = 0; i < g->n_arc; ++i) {
		asg_arc_t *p = &g->arc[i];
		if (p->del) continue;
		if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
			asg_arc_t *q;
			uint32_t u = mark[p->ul>>32^1]^1;
			int l = ug->u.a[u>>1].len - p->ol;
			if (l < 0) l = 1;
			q = asg_arc_pushp(ug->g);
			q->ol = p->ol, q->del = 0;
			q->ul = (uint64_t)u<<32 | l;
			q->v = mark[p->v];
		}
	}

    // add arcs for circular utg
    // (the same as above)
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;  // recycle `mark`
	for (i = 0; i < ug->u.n; ++i) {  // collect start/end
		if (!ug->u.a[i].circ) continue;
		mark[ug->u.a[i].start] = i<<1 | 0;
		mark[ug->u.a[i].end] = i<<1 | 1;
	}
    for (i = 0; i < g->n_arc; ++i) {  // find edges
		asg_arc_t *p = &g->arc[i];
		if (p->del) continue;
		if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
			asg_arc_t *q;
			uint32_t u = mark[p->ul>>32^1]^1;
			int l = ug->u.a[u>>1].len - p->ol;
			if (l < 0) l = 1;
			q = asg_arc_pushp(ug->g);
			q->ol = p->ol, q->del = 0;
			q->ul = (uint64_t)u<<32 | l;
			q->v = mark[p->v];
		}
	}
    
	for (i = 0; i < ug->u.n; ++i)
		asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
	asg_cleanup(ug->g);
	free(mark);
	return ug;
}

static char comp_tab[] = { // complement base
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};


void recover_fake_read(UC_Read* result, UC_Read* tmp, ma_utg_t *u,
All_reads *RNF, const ma_sub_t *coverage_cut)
{
    uint32_t j, k, ori, start, l = 0, eLen, rId, readLen;
    if(result->size < u->len)
    {
        result->size = u->len;
        result->seq = (char*)realloc(result->seq,sizeof(char)*(result->size));
    }
    char* readS = NULL;

    for (j = 0; j < u->n; ++j) {
        rId = u->a[j]>>33;
        ///uId = i;
        ori = u->a[j]>>32&1;
        start = l;
        eLen = (uint32_t)u->a[j];
        l += eLen;

        if(eLen == 0) continue;

        recover_UC_Read(tmp, RNF, rId);
        readS = tmp->seq + coverage_cut[rId].s;
        readLen = coverage_cut[rId].e - coverage_cut[rId].s;
        
        
        if (!ori) // forward strand
        {
            for (k = 0; k < eLen; k++)
            {
                result->seq[start + k] = readS[k];
            }
        }
        else
        {
            for (k = 0; k < eLen; k++)
            {
                uint8_t c = (uint8_t)readS[readLen - 1 - k];
                result->seq[start + k] = c >= 128? 'N' : comp_tab[c];
            }
        }
	}
    
    result->length = u->end - u->start;
    for (j = u->start; j < u->end; j++)
    {
        result->seq[j - u->start] = result->seq[j];
    }
}

void get_overlapLen(uint32_t rId, ma_hit_t_alloc* sources, uint32_t* exactLen, uint32_t* inexactLen)
{
    (*exactLen) = (*inexactLen) = 0;
    ma_hit_t_alloc* x = &(sources[rId]);
    ma_hit_t *h = NULL;
    uint32_t i, len;
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        len = Get_qe((*h)) - Get_qs((*h));
        if(h->el == 1)
        {
            (*exactLen) += len;
        }
        else
        {
            (*inexactLen) += len;
        }
    }
}


void reduce_ma_utg_t(ma_utg_t* collection, asg_t* read_g, ma_hit_t_alloc* sources, 
ma_sub_t *coverage_cut, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp)
{
    asg_arc_t t_f;
    asg_arc_t* av = NULL;
    uint32_t m = 0, k, i, nv;
    for (i = 0; i < collection->n; i++)
    {
        if(collection->a[i] == (uint64_t)-1) continue;
        collection->a[m] = collection->a[i];
        m++;
    }
    collection->n = m;


    uint32_t totalLen = 0, v, w, l;
    for (i = 0; i < collection->n - 1; i++)
    {
        v = (uint64_t)(collection->a[i])>>32;
        w = (uint64_t)(collection->a[i + 1])>>32; 




        /*******************************for debug************************************/
        l = (uint32_t)-1;
        av = asg_arc_a(read_g, v);
        nv = asg_arc_n(read_g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == w) 
            {
                l = asg_arc_len(av[k]);
                break;
            }
        }

        if(k == nv) 
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                {
                    l = asg_arc_len(edge->a.a[k]);
                    break;
                }
            }

            if(k == edge->a.n)
            {
                if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, v, w, &t_f)==0)
                {
                    fprintf(stderr, "####ERROR1: v>>1: %u, v&1: %u, w>>1: %u, w&1: %u, r_seq: %u\n",
                    v>>1, v&1, w>>1, w&1, read_g->r_seq);
                }
                l = asg_arc_len(t_f);
            }
        }
        if(l == (uint32_t)-1) fprintf(stderr, "[E::%s]ERROR\n", __func__);
        /*******************************for debug************************************/









        collection->a[i] = v; collection->a[i] = collection->a[i]<<32; 
        collection->a[i] = collection->a[i] | (uint64_t)(l);
        totalLen += l;
    }

    if(i < collection->n)
    {
        if(collection->circ)
        {
            v = (uint64_t)(collection->a[i])>>32;
            w = (uint64_t)(collection->a[0])>>32; 
            



            /*******************************for debug************************************/
            l = (uint32_t)-1;
            av = asg_arc_a(read_g, v);
            nv = asg_arc_n(read_g, v);
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }
            
            if(k == nv) 
            {
                for (k = 0; k < edge->a.n; k++)
                {
                    if(edge->a.a[k].del) continue;
                    if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
                    {
                        l = asg_arc_len(edge->a.a[k]);
                        break;
                    }
                }

                if(k == edge->a.n)
                {
                    if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, v, w, &t_f)==0)
                    {
                        fprintf(stderr, "####ERROR1: v>>1: %u, v&1: %u, w>>1: %u, w&1: %u, r_seq: %u\n",
                        v>>1, v&1, w>>1, w&1, read_g->r_seq);
                    }
                    l = asg_arc_len(t_f);
                }
            }
            if(l == (uint32_t)-1) fprintf(stderr, "[E::%s]ERROR\n", __func__);
            /*******************************for debug************************************/












            collection->a[i] = v; collection->a[i] = collection->a[i]<<32; 
            collection->a[i] = collection->a[i] | (uint64_t)(l);

            totalLen += l;
        }
        else
        {
            v = (uint64_t)(collection->a[i])>>32;
            l = read_g->seq[v>>1].len;
            collection->a[i] = v;
            collection->a[i] = collection->a[i]<<32;
            collection->a[i] = collection->a[i] | (uint64_t)(l);
            totalLen += l;
        }
    }

    collection->len = totalLen;
    if(!collection->circ)
    {
        collection->start = collection->a[0]>>32;
        collection->end = (collection->a[collection->n-1]>>32)^1;
    }
}

uint32_t polish_unitig_back(ma_utg_t* collection, asg_t* read_g, ma_hit_t_alloc* sources, 
ma_sub_t *coverage_cut, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp)
{
    if(collection->m == 0) return 0;
    if(collection->n < 3) return 0;
    uint32_t i, k, v, pre, afte, nv, exactLen, inexactLen, tmp_exactLen, tmp_inexactLen, skip = 0;
    asg_arc_t* av = NULL;
    asg_arc_t *pE = NULL, *aE = NULL;
    asg_arc_t t_f, t_b;
    
    for (i = 1; i < collection->n - 1; i++)
    {
        v = (uint64_t)(collection->a[i])>>32;
        pre = (uint64_t)(collection->a[i-1])>>32;
        afte = (uint64_t)(collection->a[i+1])>>32;
        if(v == (uint32_t)-1) continue;
        if(pre == (uint32_t)-1) continue;
        if(afte == (uint32_t)-1) continue;

        av = asg_arc_a(read_g, v^1);
        nv = asg_arc_n(read_g, v^1);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == (pre^1)) 
            {
                pE = &(av[k]);
                break;
            }
        }
        if(k == nv)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == (v^1) && edge->a.a[k].v == (pre^1))
                {
                    pE = &(edge->a.a[k]);
                    break;
                }
            }

            if(k == edge->a.n) fprintf(stderr, "[E::%s]ERROR 1\n", __func__);
        } 

        av = asg_arc_a(read_g, v);
        nv = asg_arc_n(read_g, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == afte)
            {
                aE = &(av[k]);
                break;
            }
        }

        if(k == nv)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == afte)
                {
                    aE = &(edge->a.a[k]);
                    break;
                }
            }
            if(k == edge->a.n) fprintf(stderr, "[E::%s]ERROR 2\n", __func__);
        }

        if(pE->el == 1 && aE->el == 1) continue;
        
        if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, pre, 
            afte, &t_f) == 0)
        {
            continue;
        }

        if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, afte^1, 
            pre^1, &t_b) == 0)
        {
            continue;
        }

        if(t_f.el == 0 || t_b.el == 0) continue;

        get_overlapLen(v>>1, sources, &exactLen, &inexactLen);
        if(pE->el == 0)
        {
            get_overlapLen(pre>>1, sources, &tmp_exactLen, &tmp_inexactLen);
            if(inexactLen < tmp_inexactLen) continue;
            if(inexactLen == tmp_inexactLen && exactLen > tmp_exactLen) continue;
        }

        if(aE->el == 0)
        {
            get_overlapLen(afte>>1, sources, &tmp_exactLen, &tmp_inexactLen);
            if(inexactLen < tmp_inexactLen) continue;
            if(inexactLen == tmp_inexactLen && exactLen > tmp_exactLen) continue;
        }

        collection->a[i] = (uint64_t)-1;
        skip++;
    }

    if(skip == 0) return 0;
    reduce_ma_utg_t(collection, read_g, sources, coverage_cut, edge, max_hang, min_ovlp);

    return 1;
}

uint32_t detect_exact_ovec(ma_utg_t* collection, asg_t* read_g, ma_hit_t_alloc* sources, 
ma_sub_t *coverage_cut, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp,
uint32_t src, uint32_t dest_idx)
{
    asg_arc_t *t = NULL, i_t;
    uint32_t i, k, dest;
    for (i = dest_idx; i < collection->n; i++)
    {
        t = NULL;
        dest = (uint64_t)(collection->a[i])>>32;
        if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, src, dest, &i_t) == 0)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == src && edge->a.a[k].v == dest)
                {
                    t = &(edge->a.a[k]);
                    break;
                }
            }
        }
        else
        {
            t = &i_t;
        }

        if(t == NULL) return (uint32_t)-1;
        if(t->el != 1) continue;
        return i;
    }

    return (uint32_t)-1;
}

int get_specific_edge(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, kvec_asg_arc_t_warp* edge, asg_t* read_g, int max_hang, 
int min_ovlp, uint32_t query, uint32_t target, asg_arc_t* t)
{
    uint32_t nv, k;
    (*t).ul = (uint64_t)-1; (*t).v = (uint32_t)-1;
    if(read_g)
    {
        asg_arc_t* av = asg_arc_a(read_g, query);
        nv = asg_arc_n(read_g, query);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == target) 
            {
                (*t) = av[k];
                break;
            }
        }
    }

    if((*t).ul == (uint64_t)-1)
    {
        if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, query, target, t)==0)
        {
            (*t).ul = (uint64_t)-1;
        }
    }
    
    if((*t).ul == (uint64_t)-1)
    {
        for (k = 0; k < edge->a.n; k++)
        {
            if(edge->a.a[k].del) continue;
            if((edge->a.a[k].ul>>32) == query && edge->a.a[k].v == target)
            {
                (*t) = edge->a.a[k];
                break;
            }
        }
        if(k == edge->a.n){
            return -1;
            //fprintf(stderr, "[E::%s] thing missing when shouldn't\n", __func__);
        } 
    }
    return 0;
}

uint32_t polish_unitig(ma_utg_t* collection, asg_t* read_g, ma_hit_t_alloc* sources, 
ma_sub_t *coverage_cut, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp)
{
    if(collection->m == 0) return 0;
    if(collection->n < 3) return 0;
    uint32_t i, k, v, pre, pre_i, afte, afte_i, exactLen, inexactLen, skip = 0;
    uint32_t min_inexactLen, max_exactLen;
    asg_arc_t pE, aE;
    pre = (uint64_t)(collection->a[0])>>32; pre_i = 0; afte_i = (uint32_t)-1;


    for (i = 1; i < collection->n - 1; i++)
    {
        if(collection->a[i] == (uint64_t)-1) continue;
        ///v and after must be available
        v = (uint64_t)(collection->a[i])>>32;
        afte = (uint64_t)(collection->a[i+1])>>32;

        int san = get_specific_edge(sources, coverage_cut, NULL, edge, pre_i == i-1? read_g:NULL, max_hang, min_ovlp, 
        v^1, pre^1, &pE);
        if (san<0) continue;
        san = get_specific_edge(sources, coverage_cut, NULL, edge, read_g, max_hang, min_ovlp, 
        v, afte, &aE);
        if (san<0) continue;
        
        if(pE.el == 1 && aE.el == 1)
        {
            pre = (uint64_t)(collection->a[i])>>32; pre_i = i;
            continue;
        } 

        ///pre must be a good read, we need to find a good after
        ///update a new afte
        afte_i = detect_exact_ovec(collection, read_g, sources, coverage_cut, edge, 
        max_hang, min_ovlp, pre, i+1);
        if(afte_i == (uint32_t)-1)
        {
            pre = (uint64_t)(collection->a[i])>>32; pre_i = i;
            continue;
        } 
        afte = (uint64_t)(collection->a[afte_i])>>32;


        min_inexactLen = (uint32_t)-1;max_exactLen = 0;
        for (k = i; k < afte_i; k++)
        {
            get_overlapLen((uint64_t)(collection->a[k])>>33, sources, &exactLen, &inexactLen);
            if(inexactLen < min_inexactLen)
            {
                min_inexactLen = inexactLen;
                max_exactLen = exactLen;
            }
        }

        get_overlapLen(pre>>1, sources, &exactLen, &inexactLen);
        if((inexactLen > min_inexactLen) || (inexactLen == min_inexactLen && exactLen <= max_exactLen))
        {
            pre = (uint64_t)(collection->a[i])>>32; pre_i = i;
            continue;
        } 

        get_overlapLen(afte>>1, sources, &exactLen, &inexactLen);
        if((inexactLen > min_inexactLen) || (inexactLen == min_inexactLen && exactLen <= max_exactLen))
        {
            pre = (uint64_t)(collection->a[i])>>32; pre_i = i;
            continue;
        } 

        for (k = i; k < afte_i; k++)
        {
            collection->a[k] = (uint64_t)-1;
            skip++;
        }
    }

    if(skip == 0) return 0;
    reduce_ma_utg_t(collection, read_g, sources, coverage_cut, edge, max_hang, min_ovlp);

    return 1;
}

inline int inter_interval(int a_s, int a_e, int b_s, int b_e, int* i_s, int* i_e)
{
    if(a_s > b_e || b_s > a_e) return 0;
    (*i_s) = MAX(a_s, b_s);
    (*i_e) = MIN(a_e, b_e);
    return 1;
}


void print_rough_inconsistent_sites(ma_utg_t* collection, uint32_t cur_i, uint32_t next_i, 
asg_t* read_g, All_reads *RNF, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
kvec_asg_arc_t_warp* edge, UC_Read* r_read, UC_Read* q_read, int max_hang, int min_ovlp,
uint32_t c_beg, uint32_t rate_thre, kvec_t_u32_warp* exact_count, kvec_t_u32_warp* total_count, 
const char* prefix, int uID, FILE* fp)
{
    uint32_t v, w, i, rate;
    int v_beg, v_end, v_sub_beg, v_sub_end, w_beg, w_end, w_sub_beg, w_sub_end, i_beg, i_end, j;
    asg_arc_t t;
    v = (uint64_t)(collection->a[cur_i])>>32;
    ///last element
    if(cur_i == collection->n-1 && next_i == collection->n)
    {
        if(!collection->circ)
        {
            next_i = (uint32_t)-1;
        }
        else
        {
            next_i = 0;
        }
    }


    if(next_i != ((uint32_t)-1))
    {
        w = (uint64_t)(collection->a[next_i])>>32;

        get_specific_edge(sources, coverage_cut, NULL, edge, read_g, max_hang, min_ovlp, v, w, &t);
        v_beg = 0; v_end = asg_arc_len(t) - 1; 
        if(v&1)
        {
            v_beg = Get_READ_LENGTH((*RNF), (v>>1)) - v_beg - 1;
            v_end = Get_READ_LENGTH((*RNF), (v>>1)) - v_end - 1;
            w = v_beg; v_beg = v_end; v_end = w;
        }
    }
    else
    {
        v_beg = 0; v_end = Get_READ_LENGTH((*RNF), (v>>1)) - 1;
    }

    kv_resize(uint32_t, exact_count->a, (uint32_t)(v_end - v_beg + 1)); 
    memset(exact_count->a.a, 0, (v_end-v_beg+1)*sizeof(uint32_t));
    kv_resize(uint32_t, total_count->a, (uint32_t)(v_end - v_beg + 1)); 
    memset(total_count->a.a, 0, (v_end-v_beg+1)*sizeof(uint32_t));
    exact_count->a.n = total_count->a.n = (v_end-v_beg+1);
    
    recover_UC_sub_Read(r_read, v_beg, v_end - v_beg + 1, 0, RNF, v>>1);
    ma_hit_t_alloc* x = &(sources[v>>1]);
    ma_hit_t *h = NULL;
    ///[v_beg, v_end] must be the end of read, which means v_beg = 0 or v_end = Get_READ_LENGTH((*RNF), (v>>1)) - 1

    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        if(inter_interval(v_beg, v_end, Get_qs((*h)), Get_qe((*h)) - 1, 
                                                &v_sub_beg, &v_sub_end) == 0)
        {
            continue;
        }

        if(h->el)
        {
            for (j = v_sub_beg; j <= v_sub_end; j++)
            {
                exact_count->a.a[j-v_beg]++;
                total_count->a.a[j-v_beg]++;
            }            
            continue;
        }

        w_beg = Get_ts((*h)); w_end = Get_te((*h)) - 1;
        if(h->rev)
        {
            w_beg = Get_READ_LENGTH((*RNF), Get_tn((*h))) - w_beg - 1;
            w_end = Get_READ_LENGTH((*RNF), Get_tn((*h))) - w_end - 1;
            w = w_beg; w_beg = w_end; w_end = w;
        }



        w_sub_beg = w_beg + (v_sub_beg - Get_qs((*h)));
        if(w_sub_beg >= (int)(Get_READ_LENGTH((*RNF), Get_tn((*h)))))
        {
            w_sub_beg = Get_READ_LENGTH((*RNF), Get_tn((*h))) - 1;
        }


        w_sub_end = w_end - ((int)(Get_qe((*h))) - 1 - v_sub_end);
        if(w_sub_end < 0) w_sub_end = 0;
        if(w_sub_beg > w_sub_end || (v_sub_end-v_sub_beg) != (w_sub_end-w_sub_beg))
        {
            for (j = v_sub_beg; j <= v_sub_end; j++)
            {
                total_count->a.a[j-v_beg]++;
            }  
            continue;
        } 

        recover_UC_sub_Read(q_read, w_sub_beg, w_sub_end-w_sub_beg +1, h->rev, RNF, Get_tn((*h)));
        if(if_exact_match(r_read->seq, r_read->length, q_read->seq, q_read->length, 
        v_sub_beg-v_beg, v_sub_end-v_beg, 0, q_read->length-1))
        {
            for (j = v_sub_beg; j <= v_sub_end; j++)
            {
                exact_count->a.a[j-v_beg]++;
                total_count->a.a[j-v_beg]++;
            }  
        }
        else
        {
            for (j = v_sub_beg; j <= v_sub_end; j++)
            {
                total_count->a.a[j-v_beg]++;
            } 
        }
    }
    
    i_beg = i_end = -1;
    v = (uint64_t)(collection->a[cur_i])>>32;
    uint32_t inexact = 0, total = 0; 
    for (i = 0; i < total_count->a.n; i++)
    {
        if(total_count->a.a[i] == 0)
        {
            rate = 100;
        }
        else
        {
            rate = ((total_count->a.a[i] - exact_count->a.a[i])*100)/total_count->a.a[i];
        }
         
        if(rate >= rate_thre)
        {
            ///start a new interval
            if(i_beg == -1 && i_end == -1)
            {
                i_beg = i_end = i;
            }
            else///extend current interval
            {
                i_end++;
            }
            total = total + total_count->a.a[i];
            inexact = inexact + (total_count->a.a[i] - exact_count->a.a[i]);
        }
        else///end an interval
        {
            if(i_beg != -1 && i_end != -1)
            {
                ///i_beg and i_end are the offsets in comparsion to v_beg
                v_sub_beg = i_beg + v_beg;
                v_sub_end = i_end + v_beg;
                if(v&1)
                {
                    i_beg = (v_end - v_beg + 1) - i_beg - 1;
                    i_end = (v_end - v_beg + 1) - i_end - 1;
                    w = i_beg; i_beg = i_end; i_end = w;
                }
                i_end++;
                rate = (total == 0)? 100 : (inexact*100)/total;

                fprintf(fp,"%s%.6d%c\t%u\t%u\t%u\t", prefix, uID, "lc"[collection->circ], 
                (uint32_t)(i_beg + c_beg), (uint32_t)(i_end + c_beg), rate);
                fprintf(fp,"%.*s", (int)Get_NAME_LENGTH((*RNF), (v>>1)), Get_NAME((*RNF), (v>>1)));
                for (j = 0; j < (int)x->length; j++)
                {
                    h = &(x->buffer[j]);
                    if(inter_interval(v_sub_beg, v_sub_end, Get_qs((*h)), Get_qe((*h)) - 1, 
                                                            &w_sub_beg, &w_sub_end) == 0)
                    {
                        continue;
                    }
                    fprintf(fp,",%.*s", (int)Get_NAME_LENGTH((*RNF), Get_tn((*h))), Get_NAME((*RNF), Get_tn((*h))));
                }
                fprintf(fp,"\n");
            }

            i_beg = i_end = -1;
            inexact = total = 0;
        }
    }

    if(i_beg != -1 && i_end != -1)
    {
        ///i_beg and i_end are the offsets in comparsion to v_beg
        v_sub_beg = i_beg + v_beg;
        v_sub_end = i_end + v_beg;
        if(v&1)
        {
            i_beg = (v_end - v_beg + 1) - i_beg - 1;
            i_end = (v_end - v_beg + 1) - i_end - 1;
            w = i_beg; i_beg = i_end; i_end = w;
        }
        i_end++;
        rate = (total == 0)? 100 : (inexact*100)/total;

        fprintf(fp,"%s%.6d%c\t%u\t%u\t%u\t", prefix, uID, "lc"[collection->circ], 
        (uint32_t)(i_beg + c_beg), (uint32_t)(i_end + c_beg), rate);
        fprintf(fp,"%.*s", (int)Get_NAME_LENGTH((*RNF), (v>>1)), Get_NAME((*RNF), (v>>1)));
        for (j = 0; j < (int)x->length; j++)
        {
            h = &(x->buffer[j]);
            if(inter_interval(v_sub_beg, v_sub_end, Get_qs((*h)), Get_qe((*h)) - 1, 
                                                    &w_sub_beg, &w_sub_end) == 0)
            {
                continue;
            }
            fprintf(fp,",%.*s", (int)Get_NAME_LENGTH((*RNF), Get_tn((*h))), Get_NAME((*RNF), Get_tn((*h))));
        }
        fprintf(fp,"\n");
    }
    
}



int get_consensus_rate(ma_utg_t* collection, uint32_t cur_i, uint32_t next_i, 
asg_t* read_g, All_reads *RNF, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
kvec_asg_arc_t_warp* edge, UC_Read* r_read, UC_Read* q_read, int max_hang, int min_ovlp,
int* r_match, int* r_total)
{
    (*r_match) = (*r_total) = 0;
    if(cur_i < 1) return -1;
    asg_arc_t *t = NULL, i_t;
    uint32_t v, w, k, v_beg, v_end, w_beg, w_end;
    int j;
    if(collection->a[cur_i] == (uint64_t)-1 || collection->a[next_i] == (uint64_t)-1) return 0;

    v = (uint64_t)(collection->a[cur_i])>>32;
    w = (uint64_t)(collection->a[next_i])>>32;

    if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, v, w, &i_t) == 0)
    {
        for (k = 0; k < edge->a.n; k++)
        {
            if(edge->a.a[k].del) continue;
            if((edge->a.a[k].ul>>32) == v && edge->a.a[k].v == w)
            {
                t = &(edge->a.a[k]);
                break;
            }
        }
    }
    else
    {
        t = &i_t;
    }

    if(t == NULL) return -1;
    v_beg = 0; v_end = asg_arc_len(*t) - 1; r_read->length = 0;
    ///cur_i must >= 1
    for (j = cur_i-1; j >= 0; j--)
    {
        if(collection->a[j] == (uint64_t)-1) continue;

        w = (uint64_t)(collection->a[j])>>32;
        t = NULL;

        if(get_edge_from_source(sources, coverage_cut, NULL, max_hang, min_ovlp, w, v, &i_t) == 0)
        {
            for (k = 0; k < edge->a.n; k++)
            {
                if(edge->a.a[k].del) continue;
                if((edge->a.a[k].ul>>32) == w && edge->a.a[k].v == v)
                {
                    t = &(edge->a.a[k]);
                    break;
                }
            }
        }
        else
        {
            t = &i_t;
        }

        if(t == NULL) break;
        
        w_beg = asg_arc_len(*t); 
        w_end = MIN(((int)(w_beg + v_end)), ((int)(read_g->seq[w>>1].len-1)));
        // if(t->el == 1)
        // {
        //     if(r_read->length == 0) recover_UC_sub_Read(r_read, v_beg, v_end - v_beg + 1, v&1, RNF, v>>1);
        //     ///for debug: check if v[v_beg, v_end] == w[w_beg, w_end]
        //     recover_UC_sub_Read(q_read, w_beg, w_end - w_beg +1, w&1, RNF, w>>1);
        //     ///q_read->length<=r_read->length
        //     if(if_exact_match(r_read->seq, r_read->length, q_read->seq, q_read->length, 
        //     0, q_read->length-1, 0, q_read->length-1) == 0)
        //     {
        //         fprintf(stderr, "ERROR: cur_i: %u, j: %d, r_len: %lld, q_len: %lld\n", cur_i, j, r_read->length, q_read->length);
        //     }
        // }
        
        //filter overlaps which cannot cover the whole [v_beg, v_end]
        if(w_end - w_beg != v_end - v_beg) break;

        (*r_total)++;

        if(t->el == 1)
        {
            (*r_match)++;
        }
        else
        {
            if(r_read->length == 0) recover_UC_sub_Read(r_read, v_beg, v_end - v_beg + 1, v&1, RNF, v>>1);
            recover_UC_sub_Read(q_read, w_beg, w_end - w_beg +1, w&1, RNF, w>>1);
            ///q_read->length<=r_read->length
            if(if_exact_match(r_read->seq, r_read->length, q_read->seq, q_read->length, 
            0, q_read->length-1, 0, q_read->length-1) == 1)
            {
                (*r_match)++;
            }
        }
    }


    return 1;
}

uint32_t polish_unitig_advance(ma_utg_t* collection, asg_t* read_g, All_reads *RNF, 
ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, kvec_asg_arc_t_warp* edge, 
UC_Read* r_read, UC_Read* q_read, int max_hang, int min_ovlp)
{
    if(collection->m == 0) return 0;
    if(collection->n < 3) return 0;
    uint32_t i, skip = 0;
    int match_v, total_v, max_i, match_max, k;
    double match_rate, match_rate_max;
        
    ///we should be able to handle i = collection->n-1
    for (i = 1; i < collection->n-1; i++)
    {
        ///in practice, collection->a[i] and collection->a[i+1] must be available
        ///collection->a[index] might be unavailable only if index < i
        if(get_consensus_rate(collection, i, i+1, read_g, RNF, sources, coverage_cut, 
        edge, r_read, q_read, max_hang, min_ovlp, &match_v, &total_v) != 1)
        {
            continue;
        }

        match_rate = (total_v == 0)? 0:((double)(match_v)/(double)(total_v));
        ///most reads support collection[i], so it is right
        if(match_v >= total_v * 0.5 && total_v > 0 && match_v > 0) continue;

        max_i = i; match_max = match_v; match_rate_max = match_rate;
        for (k = i - 1; k >= 0; k--)
        {
            if(collection->a[k] == (uint64_t)-1) continue;
            ///collection->a[k] might be unavailable, while collection->a[i+1] must be available
            ///return -1 means there is no overlap from k to i+1
            if(get_consensus_rate(collection, k, i+1, read_g, RNF, sources, coverage_cut, 
            edge, r_read, q_read, max_hang, min_ovlp, &match_v, &total_v) < 0)
            {
                break;
            }

            ///no read support k to i+1
            if(total_v == 0) break;

            match_rate = (total_v == 0)? 0:((double)(match_v)/(double)(total_v));
            if(match_rate > match_rate_max || (match_rate == match_rate_max && match_v > match_max))
            {
                max_i = k; match_max = match_v; match_rate_max = match_rate;
            }
        }
        
        ///set [max_i+1, i] to be unavailable
        for (k = max_i+1; k <= (int)i; k++)
        {
            if(collection->a[k] == (uint64_t)-1) continue;
            collection->a[k] = (uint64_t)-1;
            skip++;
        }
    }

    if(skip == 0) return 1;

    reduce_ma_utg_t(collection, read_g, sources, coverage_cut, edge, max_hang, min_ovlp);

    return 1;
}


// generate unitig sequences
int ma_ug_seq(ma_ug_t *g, asg_t *read_g, All_reads *RNF, ma_sub_t *coverage_cut,
ma_hit_t_alloc* sources, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp)
{
    UC_Read g_read;
    init_UC_Read(&g_read);
    UC_Read tmp;
    init_UC_Read(&tmp);
	///utg_intv_t *tmp;
	uint32_t i, j, k;
    uint32_t rId, /**uId,**/ori, start, eLen, readLen;
    char* readS = NULL;

    
    ///why we need n_read here? it is just beacuse one read can only be included in one untig
    ///but it is not true
	///tmp = (utg_intv_t*)calloc(n_read, sizeof(utg_intv_t));
	///number of unitigs
	for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
        if(u->m == 0) continue;
        polish_unitig(u, read_g, sources, coverage_cut, edge, max_hang, min_ovlp);
        polish_unitig_advance(u, read_g, RNF, sources, coverage_cut, edge, &g_read, &tmp, max_hang, min_ovlp);
        g->g->seq[i].len = u->len;

		uint32_t l = 0;
		u->s = (char*)calloc(1, u->len + 1);
		memset(u->s, 'N', u->len);
		for (j = 0; j < u->n; ++j) {
            rId = u->a[j]>>33;
            ///uId = i;
            ori = u->a[j]>>32&1;
            start = l;
            eLen = (uint32_t)u->a[j];
			l += eLen;

            if(eLen == 0) continue;
            if(rId < read_g->r_seq)
            {
                recover_UC_Read(&g_read, RNF, rId);
            }
            else
            {
                recover_fake_read(&g_read, &tmp, &(read_g->F_seq[rId-read_g->r_seq]),
                RNF, coverage_cut);
            }

            readS = g_read.seq + coverage_cut[rId].s;
            readLen = coverage_cut[rId].e - coverage_cut[rId].s;
            
            if (!ori) // forward strand
            {
                for (k = 0; k < eLen; k++)
                {
                    u->s[start + k] = readS[k];
                }
            }
            else
            {
                for (k = 0; k < eLen; k++)
                {
                    uint8_t c = (uint8_t)readS[readLen - 1 - k];
                    u->s[start + k] = c >= 128? 'N' : comp_tab[c];
                }
            }
		}
	}

    destory_UC_Read(&g_read);
    destory_UC_Read(&tmp);


    uint32_t n_vtx = g->g->n_seq * 2, v, nv;
    uint32_t vLen = 0;
    asg_arc_t* av = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->g->seq[v>>1].del) continue;
        av = asg_arc_a(g->g, v);
        nv = asg_arc_n(g->g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            vLen = g->g->seq[(av[i].ul>>33)].len - av[i].ol;
            av[i].ul = (av[i].ul>>32)<<32;
            av[i].ul = av[i].ul | vLen;
        }
    }

	return 0;
}


uint32_t get_ug_coverage(ma_utg_t* u, asg_t* read_g, const ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, uint8_t* r_flag)
{
    uint32_t k, j, rId, tn, is_Unitig;
    long long R_bases = 0, C_bases = 0;
    long long ret;
    ma_hit_t *h;
    if(u->m == 0) return 0;

    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        r_flag[rId] = 1;
    }

    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        R_bases += (coverage_cut[rId].e - coverage_cut[rId].s);  // if asm_opt.min_overlap_coverage is default (zero), this simply is the length of the read.
        for (j = 0; j < (uint64_t)(sources[rId].length); j++)
        {
            h = &(sources[rId].buffer[j]);
            if(h->el != 1) continue;
            tn = Get_tn((*h));
            if(read_g->seq[tn].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }
            if(read_g->seq[tn].del == 1) continue;
            if(r_flag[tn] != 1) continue;
            C_bases += (Get_qe((*h)) - Get_qs((*h)));
        }
    }


    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        r_flag[rId] = 0;
    }

    ret = C_bases/R_bases;
    ret = ret==0? 1 :ret;  // waterproof division by zero
    return ret;
}

void hamt_ug_seq_simple(asg_t *sg, ma_ug_t *ug, ma_sub_t* coverage_cut, 
                        ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp){
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    kv_destroy(new_rtg_edges.a);
}

void ma_ug_print2(const ma_ug_t *ug, All_reads *RNF, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, int print_seq, const char* prefix, FILE *fp)
{
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
        sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
		if (print_seq) fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:%u\n", name, p->s? p->s : "*", p->len,  // note: vanilla hifiasm uses rd:i instead of dp:f. The dp:f was modified for bandage.
        get_ug_coverage(p, read_g, coverage_cut, sources, ruIndex, primary_flag));
		else fprintf(fp, "S\t%s\t*\tLN:i:%d\tdp:f:%u\n", name, p->len, 
        get_ug_coverage(p, read_g, coverage_cut, sources, ruIndex, primary_flag));
        // if (print_seq) fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, p->s? p->s : "*", p->len);
		// else fprintf(fp, "S\t%s\t*\tLN:i:%d\n", name, p->len);
        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if(x<RNF->total_reads)
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
                (int)Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), 
                coverage_cut[x].s, coverage_cut[x].e, x, "apmaaa"[RNF->trio_flag[x]]);
            }
            else
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
					"FAKE", coverage_cut[x].s, coverage_cut[x].e, x, '*');
            }
            
            
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
        prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
		prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}
    #if 0
    asg_arc_t* au = NULL;
    uint32_t nu, u, v;
    for (i = 0; i < ug->u.n; ++i) {
        if(ug->u.a[i].m == 0) continue;
        if(ug->u.a[i].circ)
        {
            fprintf(fp, "L\t%s%.6dc\t+\t%s%.6dc\t+\t%dM\tL1:i:%d\n", 
            prefix, i+1, prefix, i+1, 0, ug->u.a[i].len);
            fprintf(fp, "L\t%s%.6dc\t-\t%s%.6dc\t-\t%dM\tL1:i:%d\n", 
            prefix, i+1, prefix, i+1, 0, ug->u.a[i].len);
        } 
        u = i<<1;
        au = asg_arc_a(ug->g, u);
        nu = asg_arc_n(ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], au[j].ol, asg_arc_len(au[j]));
        }


        u = (i<<1) + 1;
        au = asg_arc_a(ug->g, u);
        nu = asg_arc_n(ug->g, u);
        for (j = 0; j < nu; j++)
        {
            if(au[j].del) continue;
            v = au[j].v;
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], au[j].ol, asg_arc_len(au[j]));
        }
    }
    #endif
    free(primary_flag);    
}

void ma_ug_print(const ma_ug_t *ug, All_reads *RNF, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, const char* prefix, FILE *fp)
{
	ma_ug_print2(ug, RNF, read_g, coverage_cut, sources, ruIndex, 1, prefix, fp);
}

void hamt_ug_get_subgaware_name(char *name, int name_l, uint32_t i, const ma_ug_t *ug, const char* prefix){
    // (deprecated)
    // FUNC
    //     Incorporate trailing subgraph IDs into utg/ctg name.
    // NOTE
    //     i is the unitig index (without the direction bit)
    int verbose = 0;
    
    uint32_t readID = ug->u.a[i].a[0]>>33;
    ma_utg_subglabels_t *subgp = &R_INF.subg_label_trail->a[readID];
    ma_utg_t *p = &ug->u.a[i];
    int truncated = 0;
    
    int sancheck_label=-1, sancheck_label_length=subgp->n;

    // sancheck: all reads of a unitig should have the same subgraph label
    for (int i_read=0; i_read<ug->u.a[i].n; i_read++){
        assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].n==sancheck_label_length);  // same buffer length
        for (int i_label=0; i_label<sancheck_label_length; i_label++){
            assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].a[i_label] == subgp->a[i_label]);  // same array of labels
        }
    }

    // check if name looks too long
    if (subgp->n*3+14>name_l){
        fprintf(stderr, "[W::%s] got an extremely long trailing subg ID? will truncate", __func__);
        subgp->n = 16;
        truncated = 1;
    }
    sprintf(name, "s%d.", (int)p->subg_label);
    if (verbose) fprintf(stderr, "[debug::%s] utg %.6d, subgp n is %d\n", __func__, (int)(i)+1, (int)subgp->n);
    for (int i_sub=0; i_sub<subgp->n; i_sub++){
        sprintf(name+strlen(name), "%d", subgp->a[i_sub]);
        if (i_sub!=(subgp->n-1)) sprintf(name+strlen(name), ",");  // NOTE: name regex in gfa1.0 spec: [!-)+-<>-~][!-~]*
    }
    if (truncated){
        sprintf(name+strlen(name), "trunc");
    }
    sprintf(name+strlen(name), ".%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
}

void hamt_ug_print(const ma_ug_t *ug, All_reads *RNF, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, int print_seq, const char* prefix, FILE *fp)
{
    // with hamt extra bits
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
	uint32_t i, j, l;
	char name[32];

	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
        
        if (!asm_opt.is_use_exp_graph_cleaning){
            sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
        }else{
            sprintf(name, "s%d.%s%.6d%c", (int)p->subg_label, prefix, i + 1, "lc"[p->circ]);
        }


        // fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:%u", name, (p->s&&print_seq)? p->s : "*", p->len,  // note: vanilla hifiasm uses rd:i instead of dp:f. The dp:f was modified for bandage.
                            // get_ug_coverage(p, read_g, coverage_cut, sources, ruIndex, primary_flag));
        fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:%u", name, (p->s&&print_seq)? p->s : "*", p->len,  // note: vanilla hifiasm uses rd:i instead of dp:f. The dp:f was modified for bandage.
                            ug->utg_coverage[i]);
        if (!asm_opt.is_use_exp_graph_cleaning){
            fprintf(fp, "\n");
        }else{
            // hamt, append trailing subgraph IDs

            // sancheck: all reads of a unitig should have the same subgraph label
            uint32_t readID = ug->u.a[i].a[0]>>33;
            ma_utg_subglabels_t *subgp = &R_INF.subg_label_trail->a[readID];
            int sancheck_label=-1, sancheck_label_length=subgp->n;
            #if 0
            for (int i_read=0; i_read<ug->u.a[i].n; i_read++){
                assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].n==sancheck_label_length);  // same buffer length
                for (int i_label=0; i_label<sancheck_label_length; i_label++){
                    assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].a[i_label] == subgp->a[i_label]);  // same array of labels
                }
            }
            #endif
            
            // write
            fprintf(fp, "\tts:B:I");  // GFA1: array of uint32_t integers.
            for (int i_sub=0; i_sub<subgp->n; i_sub++){
                fprintf(fp, ",%d", (int)subgp->a[i_sub]);
            }
            fprintf(fp, "\n");
        }

        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if(x<RNF->total_reads)
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
                (int)Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), 
                coverage_cut[x].s, coverage_cut[x].e, x, "apmaaa"[RNF->trio_flag[x]]);
            }
            else
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
					"FAKE", coverage_cut[x].s, coverage_cut[x].e, x, '*');
            }
            
            
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
        if (!asm_opt.is_use_exp_graph_cleaning){
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
        }else{  // hamt
            fprintf(fp, "L\ts%d.%s%.6d%c\t%c\ts%d.%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            (int)ug->u.a[u>>1].subg_label, prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            (int)ug->u.a[v>>1].subg_label, prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
        }
	}
    free(primary_flag);    
}

int asg_cut_internal(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seq[v>>1].del) continue;
		if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
		if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_NEI) continue;
		/**
		 * so combining the last two lines, they are designed to reomve(n(1), n(2))?
							   ----->              <-------
							  |                		      |
							 n(0)--->n(1)---->n(2)---->n(3)
							  |                           |
							   ------>             <-------       
		**/
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d internal sequences\n", __func__, cnt);
	return cnt;
}





void clean_weak_ma_hit_t(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long num_sources)
{
    double startTime = Get_T();
    double tmpT;
    int placeholder1;
    long long placeholder2;
    int nb_treated = 0;

    long long i, j, index;
    uint32_t qn, tn;

    for (i = 0; i < num_sources; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);

            if(sources[i].buffer[j].del) continue;

            //if this is a weak overlap
            if(sources[i].buffer[j].ml == 0)  // is weak
            {   
                // if the current overlap of the self read is weak,
                // and we can find any strong overlap (whose target DOES have 
                // trans-overlap with the current ovlp's target)
                // for the self read, 
                // mark the current ovlp and mark-delete it in the next top code block.
                placeholder1 = check_weak_ma_hit(&(sources[qn]), reverse_sources, tn, 
                                                 Get_qs(sources[i].buffer[j]), Get_qe(sources[i].buffer[j]));
                if(!placeholder1
                /**
                ||
                !check_weak_ma_hit_reverse(&(reverse_sources[qn]), sources, tn)**/)
                {
                    sources[i].buffer[j].bl = 0;
                    index = get_specific_overlap(&(sources[tn]), tn, qn);
                    assert (index>=0);  // we have qn->tn and must be able to find its counterpart (the tn->qn overlap).
                    sources[tn].buffer[index].bl = 0;
                }
            }
        }
    }



    for (i = 0; i < num_sources; i++)
    {
        
        for (j = 0; j < sources[i].length; j++)
        {
            if(sources[i].buffer[j].del) continue;

            if(sources[i].buffer[j].bl != 0)
            {
                sources[i].buffer[j].del = 0;
            }
            else
            {
                sources[i].buffer[j].del = 1;
                nb_treated++;
            }
        }
    }

    if(VERBOSE >= 0)
    {
        fprintf(stderr, "[M::%s] treated %d, used %.1fs\n", __func__, nb_treated, Get_T()-startTime);
        fflush(stderr);
    }
}
    
typedef struct {
    long long **index;
    uint32_t **i;
    uint32_t **j;
    uint32_t *n;
    uint32_t *m;
    uint32_t i_batch_start, n_reads;
    ma_hit_t_alloc* sources; 
    ma_hit_t_alloc* reverse_sources;
}hamt_clean_weak_hit_t;

static void hamt_clean_weak_ma_hit_t_worker1(void *data, long i_r_, int tid){  // callback for kt_for()
    hamt_clean_weak_hit_t *d = (hamt_clean_weak_hit_t*)data;
    long i = i_r_ + d->i_batch_start;
    if (i>=d->n_reads) return;

    ma_hit_t_alloc* sources = d->sources; 
    ma_hit_t_alloc* reverse_sources = d->reverse_sources;

    long long j, index;
    uint32_t qn, tn, qs, qe;

    for (j = 0; j < sources[i].length; j++)
    {
        if (d->n[tid]==d->m[tid]){
            assert( (d->m[tid]>>1) < (UINT32_MAX - d->m[tid]) );  // sancheck, buffer shall not grow larger than the max possible buffer size
            d->m[tid] = d->m[tid] + (d->m[tid]>>1);
            d->index[tid] = (long long*)realloc(d->index[tid], sizeof(long long) * d->m[tid]);
            d->i[tid] = (uint32_t*)realloc(d->i[tid], sizeof(uint32_t) * d->m[tid]);
            d->j[tid] = (uint32_t*)realloc(d->j[tid], sizeof(uint32_t) * d->m[tid]);
        }

        qn = Get_qn(sources[i].buffer[j]);
        tn = Get_tn(sources[i].buffer[j]);
        if(sources[i].buffer[j].del) continue;

        if(sources[i].buffer[j].ml == 0)  // overlap is weak
        {   
            // don't do check_weak_ma_hit or it's race
            // need to log: qn, tn, qs, qe, and index of the overlap
            index = get_specific_overlap(&(sources[tn]), tn, qn);
            d->index[tid][d->n[tid]] = index;
            d->i[tid][d->n[tid]] = (uint32_t)i;
            d->j[tid][d->n[tid]] = (uint32_t)j;

            // (update linear buffer offset)
            d->n[tid] = d->n[tid]+1;
        }
    }
}
static void hamt_clean_weak_ma_hit_t_worker2(void *data, long i, int tid){  // callback for kt_for()
    hamt_clean_weak_hit_t *d = (hamt_clean_weak_hit_t*)data;
    ma_hit_t_alloc* sources = d->sources; 
    ma_hit_t_alloc* reverse_sources = d->reverse_sources;

    for (uint32_t j=0; j<sources[i].length; j++){
        if(sources[i].buffer[j].del) continue;
        sources[i].buffer[j].del = !sources[i].buffer[j].bl;
    }
}

void hamt_clean_weak_ma_hit_t(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long num_sources){
    // FUNC
    //     Threaded clean_weak_ma_hit_t
    double startTime = Get_T();
    int verbose = 0;
    int n_cpu = asm_opt.thread_num;

    hamt_clean_weak_hit_t aux;
    aux.n = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.m = (uint32_t*)calloc(n_cpu, sizeof(uint32_t));
    aux.index = (long long**)calloc(n_cpu, sizeof(long long*));
    aux.i = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));
    aux.j = (uint32_t**)calloc(n_cpu, sizeof(uint32_t*));

    aux.i_batch_start = 0;
    aux.n_reads = (uint32_t) num_sources;
    aux.sources = sources;
    aux.reverse_sources = reverse_sources;
    for (int i=0; i<n_cpu; i++){
        aux.m[i] = 1024;
        aux.index[i] = (long long*)calloc(1024, sizeof(long long));
        aux.i[i] = (uint32_t*)calloc(1024, sizeof(uint32_t));
        aux.j[i] = (uint32_t*)calloc(1024, sizeof(uint32_t));
    }

    // prepare threading in batches
    int batch_size = n_cpu * 64;
    int nb_batches = num_sources / batch_size + 1;

    double startTime_tmp;
    double accumulatedTime_a = 0;
    double accumulatedTime_b = 0;

    // step 1: collect weak hits
    uint32_t qn, tn, qs, qe, I, J;
    long long ovlp_index;
    for (int i_batch=0; i_batch<nb_batches; i_batch++){
        startTime_tmp = Get_T();
        kt_for(n_cpu, hamt_clean_weak_ma_hit_t_worker1, &aux, (long)batch_size);
        accumulatedTime_a += Get_T() - startTime_tmp;

        startTime_tmp = Get_T();
        for (int tid=0; tid<n_cpu; tid++){  // treate the sequel part
            for (int i=0; i<aux.n[tid]; i++){
                I = aux.i[tid][i];
                J = aux.j[tid][i];
                qn = Get_qn(sources[I].buffer[J]);
                tn = Get_qn(sources[I].buffer[J]);
                qs = Get_qs(sources[I].buffer[J]);
                qe = Get_qe(sources[I].buffer[J]);
                ovlp_index = aux.index[tid][i];

                if (!check_weak_ma_hit(&(sources[qn]), reverse_sources, tn, qs, qe)){
                    sources[I].buffer[J].bl = 0;
                    sources[tn].buffer[ovlp_index].bl = 0;
                }
            }
            // reset buffer
            aux.n[tid] = 0;
        }
        aux.i_batch_start += batch_size;
        accumulatedTime_b += Get_T() - startTime_tmp;
    }
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] step one finished, part1 %0.2f s, part2 %.2f s\n", __func__, 
                                           accumulatedTime_a, accumulatedTime_b);
    }
    startTime = Get_T();

    // step2: flip bits
    kt_for(n_cpu, hamt_clean_weak_ma_hit_t_worker2, &aux, (long)num_sources);

    // cleanup
    free(aux.n);
    free(aux.m);
    for (int i=0; i<n_cpu; i++) {
        free(aux.index[i]); 
        free(aux.i[i]);
        free(aux.j[i]);
    }
    free(aux.index); 
    free(aux.i); free(aux.j);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] step two finished, took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}



void debug_info_of_specfic_read(const char* name, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, int id, const char* command)
{
    long long i, j, Len;
    uint32_t tn;

    if(id == -1)
    {
        i = 0;
        Len = R_INF.total_reads;
    }
    else
    {
        i = id;
        Len = id + 1;
    }
    

    for (; i < Len; i++)
    {
        if(memcmp(name, Get_NAME(R_INF, i), Get_NAME_LENGTH(R_INF, i)) == 0)
        {
            fprintf(stderr, "\n\n\nafter %s\n", command);

            fprintf(stderr, "****************ma_hit_t (%lld)ref_read: %.*s, len: %lu****************\n", 
            i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i), Get_READ_LENGTH(R_INF, i));


            fprintf(stderr, "sources Len: %d, is_fully_corrected: %d\n", 
            sources[i].length, sources[i].is_fully_corrected);

            for (j = 0; j < sources[i].length; j++)
            {
                tn = Get_tn(sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u, ml: %u, rev: %u, el: %u, del: %u\n", 
                (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                Get_qs(sources[i].buffer[j]),
                Get_qe(sources[i].buffer[j]),
                Get_ts(sources[i].buffer[j]),
                Get_te(sources[i].buffer[j]),
                sources[i].buffer[j].ml,
                sources[i].buffer[j].rev,
                sources[i].buffer[j].el,
                (uint32_t)sources[i].buffer[j].del);
            }


            
            fprintf(stderr, "######reverse_query_read Len: %d\n", reverse_sources[i].length);



            for (j = 0; j < reverse_sources[i].length; j++)
            {
                tn = Get_tn(reverse_sources[i].buffer[j]);
                fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u, rev: %u, del: %u\n", 
                (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                Get_qs(reverse_sources[i].buffer[j]),
                Get_qe(reverse_sources[i].buffer[j]),
                Get_ts(reverse_sources[i].buffer[j]),
                Get_te(reverse_sources[i].buffer[j]),
                reverse_sources[i].buffer[j].rev,
                (uint32_t)sources[i].buffer[j].del);
            }

            break;
            
        }
    }

    fflush(stderr);

    
}

void ma_sg_print(const asg_t *g, const All_reads *RNF, const ma_sub_t *sub, FILE *fp)
{
	uint32_t i;
    for (i = 0; i < g->n_seq; ++i) 
    {
        if(!g->seq[i].del)
        {
             fprintf(fp, 
             "S\t%.*s\t*\tLN:i:%u\n",
             (int)Get_NAME_LENGTH((*RNF), i), 
             Get_NAME((*RNF), i),
             g->seq[i].len);
        }
    }

	for (i = 0; i < g->n_arc; ++i) {
		const asg_arc_t *p = &g->arc[i];
		if (sub) {
			const ma_sub_t *sq = &sub[p->ul>>33], *st = &sub[p->v>>1];

            fprintf(fp, 
            "L\t%.*s:%d-%d\t%c\t%.*s:%d-%d\t%c\t%dM\tL1:i:%u\n",   // note: `M` is `:`, changed only for bandage compat
            (int)Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            sq->s + 1, sq->e, "+-"[p->ul>>32&1],
            (int)Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            st->s + 1, st->e, "+-"[p->v&1], p->ol, (uint32_t)p->ul);
            

		} 
        else 
        {
           fprintf(fp, "L\t%.*s\t%c\t%.*s\t%c\t%dM\tL1:i:%u\n",    // note: `M` is `:`, changed only for bandage compat
            (int)Get_NAME_LENGTH((*RNF), p->ul>>33), 
            Get_NAME((*RNF), p->ul>>33),
            "+-"[p->ul>>32&1],
            (int)Get_NAME_LENGTH((*RNF), p->v>>1), 
            Get_NAME((*RNF), p->v>>1),
            "+-"[p->v&1], p->ol, (uint32_t)p->ul);
		}
	}
}


void ma_ug_print_simple(const ma_ug_t *ug, All_reads *RNF, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, const char* prefix, FILE *fp)
{
	ma_ug_print2(ug, RNF, read_g, coverage_cut, sources, ruIndex, 0, prefix, fp);
}

void ma_ug_print_bed(const ma_ug_t *g, asg_t *read_g, All_reads *RNF, ma_sub_t *coverage_cut,
ma_hit_t_alloc* sources, kvec_asg_arc_t_warp* edge, int max_hang, int min_ovlp, uint32_t rate_thres,
const char* prefix, FILE *fp)
{
	UC_Read g_read;
    init_UC_Read(&g_read);
    UC_Read tmp;
    init_UC_Read(&tmp);
    kvec_t_u32_warp exact_count, total_count;
    kv_init(exact_count.a);
    kv_init(total_count.a);
    uint32_t i, j, l, eLen, start;
    for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
        if(u->m == 0) continue;
        if(u->n < 2) continue;
        l = 0;
        for (j = 0; j < u->n; ++j) 
        {
            start = l;
            eLen = (uint32_t)u->a[j];
			l += eLen;

            print_rough_inconsistent_sites(u, j, j+1, read_g, RNF, sources, coverage_cut, 
            edge, &g_read, &tmp, max_hang, min_ovlp, start, rate_thres, &exact_count, 
            &total_count, prefix, i+1, fp);
        }
    }

    destory_UC_Read(&g_read);
    destory_UC_Read(&tmp);
    kv_destroy(exact_count.a);
    kv_destroy(total_count.a);
}


int asg_arc_cut_long_tip_primary(asg_t *g, ma_ug_t *ug, float drop_ratio)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, v_maxLen_i = (uint32_t)-1;
    long long ll, v_maxLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), flag;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        v_maxLen = -1;
        v_maxLen_i = (uint32_t)-1;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(ug == NULL)
                {
                    detect_single_path_with_dels(g, av[i].v, &convex, &ll, NULL);
                }
                else
                {
                    untig_detect_single_path_with_dels(g, ug, av[i].v, &convex, &ll, NULL);
                }
                
                
                if(v_maxLen < ll)
                {
                    v_maxLen = ll;
                    v_maxLen_i = i;
                }
            }
        }

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_maxLen_i == i) continue;

                b.b.n = 0;

                if(ug == NULL)
                {
                    flag = detect_single_path_with_dels(g, av[i].v, &convex, &ll, &b);
                }
                else
                {
                    flag = untig_detect_single_path_with_dels(g, ug, av[i].v, &convex, &ll, &b);
                }
                
                if(flag == END_TIPS)
                {
                    if(v_maxLen*drop_ratio > ll)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }

                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }
                    }
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}


int asg_arc_cut_long_tip_primary_complex(asg_t *g, float drop_ratio, uint32_t stops_threshold)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;      
        flag = detect_single_path_with_dels(g, v^1, &convex, &ll, NULL);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        convex = convex^1;ll--;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);


        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                ///if stops_threshold = 1, 
                ///detect_single_path_with_dels_n_stops() is detect_single_path_with_dels()
                detect_single_path_with_dels_n_stops(g, a_convex[i].v, &convex, &ll, &max_stopLen,
                NULL, stops_threshold);

                if(convex == v) continue;

                if(ll*drop_ratio > convexLen && max_stopLen*2>ll)
                {

                    b.b.n = 0; 
                    flag = detect_single_path_with_dels(g, v^1, &convex, &ll, &b);
                    if(b.b.n < 2) break;
                    b.b.n--;


                    n_reduced++;
                    uint64_t k;

                    for (k = 0; k < b.b.n; k++)
                    {
                        g->seq[b.b.a[k]].c = ALTER_LABLE;
                    }

                    for (k = 0; k < b.b.n; k++)
                    {
                        asg_seq_drop(g, b.b.a[k]);
                    }
                    break;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}


int untig_asg_arc_cut_long_tip_primary_complex(ma_ug_t *ug, float drop_ratio, 
uint32_t stops_threshold)
{
    double startTime = Get_T();
    asg_t *g = ug->g;
    
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t buf;
    memset(&buf, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        /****************************may have bugs********************************/
        buf.b.n = 0;       
        flag = untig_detect_single_path_with_dels(g, ug, v^1, &convex, &ll, &buf);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        get_real_length(g, convex, &convex);
        /****************************may have bugs********************************/
        convex = convex^1;
        ///note the convexLen here
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);


        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                ///if stops_threshold = 1, 
                ///untig_detect_single_path_with_dels_n_stops() is untig_detect_single_path_with_dels()
                untig_detect_single_path_with_dels_n_stops(g, ug, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);
                if(convex == v) continue;

                if(ll*drop_ratio > convexLen && max_stopLen*2>ll)
                {
                    n_reduced++;
                    uint64_t k;

                    for (k = 0; k < buf.b.n; k++)
                    {
                        g->seq[buf.b.a[k]].c = ALTER_LABLE;
                    }

                    for (k = 0; k < buf.b.n; k++)
                    {
                        asg_seq_drop(g, buf.b.a[k]);
                    }
                    break;
                }
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(buf.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}



int untig_asg_arc_cut_long_equal_tips_assembly(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i;
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), n_tips;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;
        is_hap = 0;

        ///there must be more than 1 out-edges
        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

                if(flag == END_TIPS)
                {
                    n_tips++;
                }
            }
        }

        ///at least one tip
        if(n_tips > 0)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    if(untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b) 
                    != END_TIPS)
                    {
                        continue;
                    }
                    //we can only cut tips
                    if(check_if_diploid_untigs(g, read_sg, av[base_maxLen_i].v, av[i].v, &(ug->u),
                    reverse_sources, miniedgeLen, 0.3, &b_0, &b_1, 1, ruIndex) == 1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }
        }

        if(is_hap > 0)
        {
            i = base_maxLen_i;
            b.b.n = 0;
            untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
            uint64_t k;
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]].c = HAP_LABLE;
            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

int asg_arc_cut_long_equal_tips_assembly(asg_t *g, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v), n_tips;
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

                if(flag == END_TIPS)
                {
                    n_tips++;
                }
            }
        }

        ///at least one tip
        if(n_tips > 0)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    if(detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b) 
                    != END_TIPS)
                    {
                        continue;
                    }
                    //we can only cut tips
                    /****************************may have bugs********************************/
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen, ruIndex)==1)
                    {/****************************may have bugs********************************/
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }
                        is_hap++;
                    }
                }
            }
        }


        if(is_hap > 0)
        {
            i = base_maxLen_i;
            b.b.n = 0;
            detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);

            uint64_t k;
            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]].c = HAP_LABLE;
            }
        }


    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}


int asg_arc_simple_large_bubbles(asg_t *g, ma_hit_t_alloc* reverse_sources, long long miniedgeLen, 
R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i, all_covex;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        all_covex = -1;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                if(flag != TWO_INPUT && flag != MUL_INPUT)
                {
                    break;
                }

                if(all_covex != -1 && (uint32_t)all_covex != convex)
                {
                    break;
                }

                if(all_covex == -1)
                {
                    all_covex = convex;
                }
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

            }
        }


        if(i == nv)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                    if(b.b.n < 2) continue;
                    b.b.n--;

                    //we can only cut tips
                    /****************************may have bugs********************************/
                    if(check_if_diploid(av[base_maxLen_i].v, av[i].v, g, 
                    reverse_sources, miniedgeLen, ruIndex)==1)
                    {/****************************may have bugs********************************/
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }

            if(is_hap > 0)
            {
                i = base_maxLen_i;
                b.b.n = 0;
                detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                uint64_t k;
                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]].c = HAP_LABLE;
                }
            }
        }

    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

int untig_asg_arc_simple_large_bubbles(ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap;
    long long ll, base_maxLen, base_maxLen_i, all_covex;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        all_covex = -1;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                flag = untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, NULL);
                if(flag != TWO_INPUT && flag != MUL_INPUT)
                {
                    break;
                }

                get_real_length(g, convex, &convex);

                if(all_covex != -1 && (uint32_t)all_covex != convex)
                {
                    break;
                }

                if(all_covex == -1)
                {
                    all_covex = convex;
                }
                
                if(base_maxLen < ll)
                {
                    base_maxLen = ll;
                    base_maxLen_i = i;
                }

            }
        }


        if(i == nv)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_maxLen_i) continue;
                if (!av[i].del)
                {
                    b.b.n = 0;
                    untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                    

                    //we can only cut tips
                    if(check_if_diploid_untigs(g, read_sg, av[base_maxLen_i].v, av[i].v, &(ug->u),
                    reverse_sources, miniedgeLen, 0.3, &b_0, &b_1, 1, ruIndex) == 1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }

                        is_hap++;
                    }
                }
            }

            if(is_hap > 0)
            {
                i = base_maxLen_i;
                b.b.n = 0;
                untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                uint64_t k;
                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]].c = HAP_LABLE;
                }
            }
        }

    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

uint32_t get_num_trio_flag(ma_ug_t *ug, uint32_t v, uint32_t flag)
{
    if(flag == (uint32_t)-1) return 0;
    ma_utg_t* u = NULL;
    asg_t* nsg = ug->g;
    uint32_t k, rId, flag_occ = 0;
    if (nsg->seq[v].del) return 0;
    u = &(ug->u.a[v]);
    if(u->m == 0) return 0;
    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        if(R_INF.trio_flag[rId] == flag) flag_occ++;
    }

    return flag_occ;
}

int untig_asg_arc_simple_large_bubbles_trio(ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, R_to_U* ruIndex, uint32_t positive_flag, uint32_t negative_flag)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap, k, to_replace;
    long long ll, base_maxLen, base_maxPositive, base_minNegative, base_minNonPositive, base_best_i, all_covex, curPositive, curNonPositive, curNegative;
    long long tmp, max_stop_nodeLen, max_stop_baseLen, cur_weight = 0, max_weight = 0;
    uint32_t non_positive_flag = (uint32_t)-1;
    if(positive_flag == FATHER) non_positive_flag = MOTHER;
    if(positive_flag == MOTHER) non_positive_flag = FATHER;

    buf_t buffer;
    memset(&buffer, 0, sizeof(buf_t));
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_minNegative = INT64_MAX;
        base_minNonPositive = INT64_MAX;
        base_maxPositive = -1;
        base_maxLen = -1;
        base_best_i = -1;
        all_covex = -1;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                buffer.b.n = 0;
                ///flag = untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                flag = get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, &max_stop_baseLen, 1, &buffer);
                if(flag != MUL_INPUT) break;
                // if(flag != TWO_INPUT && flag != MUL_INPUT)
                // {
                //     break;
                // }

                get_real_length(g, convex, &convex);

                if(all_covex != -1 && (uint32_t)all_covex != convex)
                {
                    break;
                }

                if(all_covex == -1)
                {
                    all_covex = convex;
                }
                

                curPositive = curNegative = curNonPositive = 0;
                for (k = 0; k < buffer.b.n; k++)
                {
                    curPositive = curPositive + get_num_trio_flag(ug, buffer.b.a[k]>>1, positive_flag);
                    curNegative = curNegative + get_num_trio_flag(ug, buffer.b.a[k]>>1, negative_flag);
                    if(non_positive_flag != (uint32_t)-1)
                    {
                        curNonPositive = curNonPositive + get_num_trio_flag(ug, buffer.b.a[k]>>1, non_positive_flag);
                    }
                }

                to_replace = 0;
                cur_weight = (long long)curPositive - ((long long)curNegative + (long long)curNonPositive);
                max_weight = (long long)base_maxPositive - ((long long)base_minNegative + (long long)base_minNonPositive);
                if(cur_weight > max_weight)
                {
                    to_replace = 1;
                }
                else if (cur_weight == max_weight)
                {
                    if(ll > base_maxLen)
                    {
                        to_replace = 1;
                    }
                }


                if(base_best_i == -1 || to_replace == 1)
                {
                    base_minNegative = curNegative;
                    base_minNonPositive = curNonPositive;
                    base_maxPositive = curPositive;
                    base_maxLen = ll;
                    base_best_i = i;
                }
            }
        }


        if(i == nv)
        {
            for (i = 0; i < nv; i++)
            {
                if(i == base_best_i) continue;
                if (!av[i].del)
                {
                    buffer.b.n = 0;
                    ///untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                    get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, &max_stop_baseLen, 1, &buffer);

                    //we can only cut tips
                    /**if(check_if_diploid_untigs(g, read_sg, av[base_best_i].v, av[i].v, &(ug->u),
                    reverse_sources, miniedgeLen, 0.3, &b_0, &b_1, 1, ruIndex) == 1)**/
                    if(check_different_haps(g, ug, read_sg, av[base_best_i].v, av[i].v, reverse_sources, &b_0, &b_1,
                    ruIndex, miniedgeLen, 1)==PLOID)
                    {
                        n_reduced++;

                        for (k = 0; k < buffer.b.n; k++)
                        {
                            g->seq[buffer.b.a[k]>>1].c = ALTER_LABLE;
                        }
                        for (k = 0; k < buffer.b.n; k++)
                        {
                            asg_seq_drop(g, buffer.b.a[k]>>1);
                        }

                        is_hap++;
                    }
                }
            }

            if(is_hap > 0)
            {
                i = base_best_i;
                buffer.b.n = 0;
                ///untig_detect_single_path_with_dels_contigLen(g, av[i].v, &convex, &ll, &b);
                get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, &max_stop_baseLen, 1, &buffer);
                for (k = 0; k < buffer.b.n; k++)
                {
                    g->seq[buffer.b.a[k]>>1].c = HAP_LABLE;
                }
            }
        }

    }


    asg_cleanup(g);
    asg_symm(g);
    free(buffer.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}



int asg_arc_cut_long_equal_tips_assembly_complex(asg_t *g, ma_hit_t_alloc* reverse_sources, 
long long miniedgeLen, uint32_t stops_threshold, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {

        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        flag = detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, NULL);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        convex = convex^1;
        ///uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = (uint32_t)-1, convex_i = (uint32_t)-1;
        asg_arc_t *a_convex = asg_arc_a(g, convex);
        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                detect_single_path_with_dels_contigLen(g, a_convex[i].v, &convex, &ll, NULL);
                if(convex == v)
                {
                    convexLen = ll;
                    convex_i = i;
                    break;
                }
            }
        }

        
        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                if(i == convex_i) continue;
                
                detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);
                ///threshold = 0.8
                if(ll > convexLen  && max_stopLen*1.25>ll)
                {
                    ///keep all nodes of this tip
                    b.b.n = 0; 
                    detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, &b);
                    if(b.b.n < 2) break;
                    b.b.n--;
                    ///keep all nodes of this tip

                     //we can only cut tips
                    if(check_if_diploid_primary_complex(v^1, a_convex[i].v, g, 
                    reverse_sources, miniedgeLen, stops_threshold, 1, ruIndex)==1)
                    {
                        n_reduced++;
                        uint64_t k;

                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < b.b.n; k++)
                        {
                            asg_seq_drop(g, b.b.a[k]);
                        }


                        ///lable the primary one
                        b.b.n = 0; 
                        detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, 
                        &ll, &max_stopLen, &b, stops_threshold);
                        for (k = 0; k < b.b.n; k++)
                        {
                            g->seq[b.b.a[k]].c = HAP_LABLE;
                        }

                        break;
                    }

                }
                
            }
        }    
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}



int untig_asg_arc_cut_long_equal_tips_assembly_complex(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll, max_stopLen;

    buf_t buf, b_0, b_1;
    memset(&buf, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;
        /****************************may have bugs********************************/
        buf.b.n = 0; 
        flag = untig_detect_single_path_with_dels_contigLen(g, v^1, &convex, &ll, &buf);
        if(flag != TWO_INPUT && flag != MUL_INPUT) continue;
        get_real_length(g, convex, &convex);
        /****************************may have bugs********************************/
        convex = convex^1;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);

        for (i = 0; i < n_convex; i++)
        {
            if (!a_convex[i].del)
            {
                untig_detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, &convex, &ll, 
                &max_stopLen, NULL, stops_threshold);

                if(convex == v) continue;

                ///threshold = 0.8
                if(ll > convexLen  && max_stopLen*1.25>ll)
                {
                    
                     //we can only cut tips
                    if(check_if_diploid_untigs_complex(g, read_sg, v^1, a_convex[i].v, 
                    &(ug->u), reverse_sources, miniedgeLen, stops_threshold,
                    0.3, &b_0, &b_1, 1, ruIndex)==1)
                    {
                        n_reduced++;
                        uint64_t k;
                        for (k = 0; k < buf.b.n; k++)
                        {
                            g->seq[buf.b.a[k]].c = ALTER_LABLE;
                        }
                        for (k = 0; k < buf.b.n; k++)
                        {
                            asg_seq_drop(g, buf.b.a[k]);
                        }




                        ///lable the primary one
                        b_0.b.n = 0; 
                        untig_detect_single_path_with_dels_contigLen_complex(g, a_convex[i].v, 
                        &convex, &ll, &max_stopLen, &b_0, stops_threshold);
                        for (k = 0; k < b_0.b.n; k++)
                        {
                            g->seq[b_0.b.a[k]].c = HAP_LABLE;
                        }
                        
                        break;
                    }

                }
                
            }
        }    
    }


    asg_cleanup(g);
    asg_symm(g);
    free(buf.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}



int untig_asg_arc_cut_chimeric_back(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, float drop_rate,
R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, w1, w2, wv, nw, n_vtx = g->n_seq * 2, n_reduced = 0, convex, convex_T;
    asg_arc_t *aw;
    long long ll, max_stopLen;

    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v = 0; v < n_vtx; ++v) 
    {

        uint32_t i;
        ///some node could be deleted
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        if(get_real_length(g, v, NULL) != 1 || get_real_length(g, v^1, NULL) != 1 ) continue;

        get_real_length(g, v, &w1);
        if(get_real_length(g, w1^1, NULL)<=1) continue;

        get_real_length(g, v^1, &w2);
        if(get_real_length(g, w2^1, NULL)<=1) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w1, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) continue;
        if(ll*0.4 > max_stopLen) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w2, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) continue;
        if(ll*0.4 > max_stopLen) continue;

        stops_threshold++;


        w1 = w1^1;wv=v^1;aw = asg_arc_a(g, w1); nw = asg_arc_n(g, w1);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1)==(v>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v>>1)) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(convex == convex_T) break;

            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        w2 = w2^1;wv=v;aw = asg_arc_a(g, w2); nw = asg_arc_n(g, w2);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < EvaluateLen((*ug).u, v>>1)) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1) == (v>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v>>1)) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(convex == convex_T) break;

            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        n_reduced++;g->seq[v>>1].c = ALTER_LABLE;asg_seq_drop(g, v>>1);
        ///fprintf(stderr, "v>>1: %u\n", v>>1);
    }


    asg_cleanup(g);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}


int untig_asg_arc_cut_chimeric(ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, float drop_rate,
R_to_U* ruIndex)
{
    asg_t *g = ug->g;
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t i, k, v_i, v_beg, v_end, selfLen, w1, w2, wv, nw, n_vtx = g->n_seq * 2, n_reduced = 0, convex, convex_T;
    asg_arc_t *aw;
    long long ll, max_stopLen;

    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v_i = 0; v_i < n_vtx; ++v_i) 
    {

        v_beg = v_i;
        ///some node could be deleted
        if (g->seq[v_beg>>1].del || g->seq[v_beg>>1].c == ALTER_LABLE) continue;
        if(get_real_length(g, v_beg, NULL) != 1) continue;

        get_real_length(g, v_beg, &w1);
        if(get_real_length(g, w1^1, NULL)<=1) continue;

        
        selfLen = get_long_tip_length(g, &(ug->u), v_beg^1, &v_end, NULL);
        if(ug->u.a[v_beg>>1].circ) continue;
        if(selfLen == (uint32_t)-1) continue;
        if(get_real_length(g, v_end, NULL) != 1) continue;


        get_real_length(g, v_end, &w2);
        if(get_real_length(g, w2^1, NULL)<=1) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w1, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < selfLen) continue;
        if(ll*0.4 > max_stopLen) continue;


        untig_detect_single_path_with_dels_n_stops(g, ug, w2, &convex, &ll, &max_stopLen, NULL, 
        stops_threshold);
        if(ll*drop_rate < selfLen) continue;
        if(ll*0.4 > max_stopLen) continue;

        stops_threshold++;


        w1 = w1^1;wv=v_beg^1;aw = asg_arc_a(g, w1); nw = asg_arc_n(g, w1);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < selfLen) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1)==(v_beg>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v_beg>>1)) continue;

            b_0.b.n = 0;
            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, &b_0, stops_threshold);
            ///if(convex == convex_T) break;
            for (k = 0; k < b_0.b.n; k++)
            {
                if(b_0.b.a[k] == (convex_T>>1)) break;
            }
            if(k != b_0.b.n) break;
            

            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        w2 = w2^1;wv=v_end^1;aw = asg_arc_a(g, w2); nw = asg_arc_n(g, w2);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, NULL, stops_threshold);
            if(ll*drop_rate < selfLen) break;
            if(ll*0.4 > max_stopLen) continue;

            if((aw[i].v>>1) == (v_end>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v_end>>1)) continue;

            b_0.b.n = 0;
            untig_detect_single_path_with_dels_n_stops(g, ug, aw[i].v, &convex, &ll, 
            &max_stopLen, &b_0, stops_threshold);
            ///if(convex == convex_T) break;
            for (k = 0; k < b_0.b.n; k++)
            {
                if(b_0.b.a[k] == (convex_T>>1)) break;
            }
            if(k != b_0.b.n) break;


            if(check_if_diploid_untigs_complex(g, read_sg, wv, aw[i].v, &(ug->u), reverse_sources, 
            miniedgeLen, stops_threshold, 0.3, &b_0, &b_1, 1, ruIndex)==1)
            {
                break;
            }
        }
        if(i!=nw) continue;






        n_reduced++;
        b_0.b.n = 0;
        get_long_tip_length(g, &(ug->u), v_beg^1, &v_end, &b_0);
        for (k = 0; k < b_0.b.n; k++)
        {
            g->seq[b_0.b.a[k]>>1].c = ALTER_LABLE;
        }

        for (k = 0; k < b_0.b.n; k++)
        {
            asg_seq_drop(g, b_0.b.a[k]>>1);
        }
    }


    asg_cleanup(g);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}


void output_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);
    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "Writing raw unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.r_utg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);
    sprintf(gfa_name, "%s.r_utg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);
    if(asm_opt.bed_inconsist_rate != 0)
    {
        sprintf(gfa_name, "%s.r_utg.lowQ.bed", output_file_name);
        output_file = fopen(gfa_name, "w");
        ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
        max_hang, min_ovlp, asm_opt.bed_inconsist_rate, "utg", output_file);
        fclose(output_file);
    }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}

void hamt_output_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, const char *suffix, const char *seq_prefix,
                              ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp, int label)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    if (label>=0){
        ug = ma_ug_gen_primary(sg, (uint8_t)label);
    }else{
        ug = ma_ug_gen(sg);
    }
    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "[M::%s] Writing GFA... \n", __func__);
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s.gfa", output_file_name, suffix);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, seq_prefix, output_file);
    fclose(output_file);
    sprintf(gfa_name, "%s.%s.noseq.gfa", output_file_name, suffix);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, seq_prefix, output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}


void hamt_output_unitig_graph_advance(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, const char *suffix, const char *seq_prefix,
                              ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp, int label)
{
    // NOTE
    //     include subgraph ID in unitig/contig name
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, label);
    hamt_ug_util_BFS_markSubgraph(ug, -1);// collect subgraph info

    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "[M::%s] Writing GFA... \n", __func__);
    
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s.gfa", output_file_name, suffix);
    FILE* output_file = fopen(gfa_name, "w");
    hamt_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, 1, seq_prefix, output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.%s.noseq.gfa", output_file_name, suffix);
    output_file = fopen(gfa_name, "w");
    hamt_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, 0, seq_prefix, output_file);
    fclose(output_file);

    free(gfa_name);
    hamt_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}

void hamt_ug_print_debug(const ma_ug_t *ug, const ma_ug_t *ug0, All_reads *RNF, asg_t* read_g, const ma_sub_t *coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, int print_seq, const char* prefix, FILE *fp)
{
    // Monkey patch, use ug->utg_coverage[i] rather than calculating the coverage on the fly.
    //  This was used to easily enable one integer labeling to visualize in bandage, no other usage.
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
	uint32_t i, j, l;
	char name[32];

	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
        
        if (!asm_opt.is_use_exp_graph_cleaning){
            sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
        }else{
            sprintf(name, "s%d.%s%.6d%c", (int)p->subg_label, prefix, i + 1, "lc"[p->circ]);
        }


        // fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:%u", name, (p->s&&print_seq)? p->s : "*", p->len,  // note: vanilla hifiasm uses rd:i instead of dp:f. The dp:f was modified for bandage.
                            // get_ug_coverage(p, read_g, coverage_cut, sources, ruIndex, primary_flag));
        fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:%u", name, (p->s&&print_seq)? p->s : "*", p->len,  // note: vanilla hifiasm uses rd:i instead of dp:f. The dp:f was modified for bandage.
                            ug0->utg_coverage[i]);
        if (!asm_opt.is_use_exp_graph_cleaning){
            fprintf(fp, "\n");
        }else{
            // hamt, append trailing subgraph IDs
            uint32_t readID = ug->u.a[i].a[0]>>33;
            ma_utg_subglabels_t *subgp = &R_INF.subg_label_trail->a[readID];
            int sancheck_label=-1, sancheck_label_length=subgp->n;
            #if 0
            // sancheck: all reads of a unitig should have the same subgraph label
            for (int i_read=0; i_read<ug->u.a[i].n; i_read++){
                assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].n==sancheck_label_length);  // same buffer length
                for (int i_label=0; i_label<sancheck_label_length; i_label++){
                    assert(R_INF.subg_label_trail->a[ug->u.a[i].a[i_read]>>33].a[i_label] == subgp->a[i_label]);  // same array of labels
                }
            }
            #endif
            
            // write
            fprintf(fp, "\tts:B:I");  // GFA1: array of uint32_t integers.
            for (int i_sub=0; i_sub<subgp->n; i_sub++){
                fprintf(fp, ",%d", (int)subgp->a[i_sub]);
            }
            fprintf(fp, "\n");
        }

        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if(x<RNF->total_reads)
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
                (int)Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), 
                coverage_cut[x].s, coverage_cut[x].e, x, "apmaaa"[RNF->trio_flag[x]]);
            }
            else
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
					"FAKE", coverage_cut[x].s, coverage_cut[x].e, x, '*');
            }
            
            
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
        if (!asm_opt.is_use_exp_graph_cleaning){
            fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
        }else{  // hamt
            fprintf(fp, "L\ts%d.%s%.6d%c\t%c\ts%d.%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
            (int)ug->u.a[u>>1].subg_label, prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
            (int)ug->u.a[v>>1].subg_label, prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
        }
	}
    free(primary_flag);    
}
void hamt_output_unitig_graph_advance_debug(asg_t *sg, ma_ug_t *ug0, ma_sub_t* coverage_cut, char* output_file_name, const char *suffix, const char *seq_prefix,
                              ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp, int label)
{
    // (companion monkey patch of the function above.)
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, label);
    hamt_ug_util_BFS_markSubgraph(ug, label);// collect subgraph info

    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "[M::%s] Writing GFA... \n", __func__);
    
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s.gfa", output_file_name, suffix);
    FILE* output_file = fopen(gfa_name, "w");
    hamt_ug_print_debug(ug, ug0, &R_INF, sg, coverage_cut, sources, ruIndex, 1, seq_prefix, output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.%s.noseq.gfa", output_file_name, suffix);
    output_file = fopen(gfa_name, "w");
    hamt_ug_print_debug(ug, ug0, &R_INF, sg, coverage_cut, sources, ruIndex, 0, seq_prefix, output_file);
    fclose(output_file);

    free(gfa_name);
    hamt_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}

void merge_unitig_content(ma_utg_t* collection, ma_ug_t* ug, asg_t* read_g, kvec_asg_arc_t_warp* edge)
{
    if(collection->m == 0) return;
    uint32_t i, j, index = 0, uId, ori;
    uint64_t totalLen;
    ma_utg_t* query = NULL;
    for (i = 0; i < collection->n; i++)
    {
        uId = collection->a[i]>>33;
        ori = collection->a[i]>>32&1;
        query = &(ug->u.a[uId]);
        index += query->n;
    }

    uint64_t* buffer = (uint64_t*)malloc(sizeof(uint64_t)*index);
    uint64_t* aim = NULL;
    index = 0;
    for (i = 0; i < collection->n; i++)
    {
        uId = collection->a[i]>>33;
        ori = collection->a[i]>>32&1;
        query = &(ug->u.a[uId]);
        aim = buffer + index;
        if(ori == 1)
        {
            for (j = 0; j < query->n; j++)
            {
                aim[query->n - j - 1] = (query->a[j])^(uint64_t)(0x100000000);
            }
        }
        else
        {
            for (j = 0; j < query->n; j++)
            {
                aim[j] = query->a[j];
            }
        }
        index += query->n;
        free(query->a);query->m=0;query->a=NULL;
    }

    if(index == 0) return;

    fill_unitig(buffer, index, read_g, edge, collection->circ, &totalLen);

    ///important. must be here
    if(collection->n == 1 && ug->u.a[collection->a[0]>>33].circ) collection->circ = 1;
    
    free(collection->a);
    collection->a = buffer;
    collection->n = collection->m = index;
    collection->len = totalLen;
    if(!collection->circ)
    {
        collection->start = collection->a[0]>>32;
        collection->end = (collection->a[collection->n-1]>>32)^1;
    }
    else
    {
        collection->start = collection->end = UINT32_MAX;
    }
    
    
}

void print_unitig(ma_utg_t* collection, ma_ug_t* ug)
{
    if(collection->m == 0) return;
    uint32_t i, index = 0, uId, ori, l;
    ma_utg_t* query = NULL;
    for (i = 0; i < collection->n; i++)
    {
        uId = collection->a[i]>>33;
        ori = collection->a[i]>>32&1;
        l = (uint32_t)(collection->a[i]);
        fprintf(stderr, "i: %u, uId: %u, ori: %u, l: %u\n", i, uId, ori, l);
        if(ug)
        {
            query = &(ug->u.a[uId]);
            index += query->n;
        }
    }
    fprintf(stderr, "index: %u\n", index);
}


void print_specfic_read_ovlp(uint32_t Rid,  ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
const char* command)
{
    long long j, i = Rid;
    uint32_t tn;


    fprintf(stderr, "\n\n\nafter %s\n", command);

    fprintf(stderr, "****************ma_hit_t (%lld)ref_read: %.*s****************\n", 
    i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i));


    fprintf(stderr, "sources Len: %d, is_fully_corrected: %d\n", 
    sources[i].length, sources[i].is_fully_corrected);

    for (j = 0; j < sources[i].length; j++)
    {
        tn = Get_tn(sources[i].buffer[j]);
        fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u, ml: %u, rev: %u, el: %u, del: %u\n", 
        (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
        Get_qs(sources[i].buffer[j]),
        Get_qe(sources[i].buffer[j]),
        Get_ts(sources[i].buffer[j]),
        Get_te(sources[i].buffer[j]),
        sources[i].buffer[j].ml,
        sources[i].buffer[j].rev,
        sources[i].buffer[j].el,
        (uint32_t)sources[i].buffer[j].del);
    }


    
    fprintf(stderr, "######reverse_query_read Len: %d\n", reverse_sources[i].length);



    for (j = 0; j < reverse_sources[i].length; j++)
    {
        tn = Get_tn(reverse_sources[i].buffer[j]);
        fprintf(stderr, "target: %.*s, qs: %u, qe: %u, ts: %u, te: %u, del: %u\n", 
        (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
        Get_qs(reverse_sources[i].buffer[j]),
        Get_qe(reverse_sources[i].buffer[j]),
        Get_ts(reverse_sources[i].buffer[j]),
        Get_te(reverse_sources[i].buffer[j]),
        (uint32_t)sources[i].buffer[j].del);
    }

    
        

    fflush(stderr);

    
}


uint32_t print_debug_gfa(asg_t *read_g, ma_ug_t *ug, ma_sub_t* coverage_cut, const char* output_file_name, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    if(ug == NULL) ug = ma_ug_gen(read_g);

    uint32_t i;
    for (i = 0; i < ug->u.n; ++i) 
    {
        ma_utg_t *u = &ug->u.a[i];
        if(u->m == 0 || ug->g->seq[i].c == ALTER_LABLE)
        {
            asg_seq_del(ug->g, i);
            if(ug->u.a[i].m!=0)
            {
                ug->u.a[i].m = ug->u.a[i].n = 0;
                free(ug->u.a[i].a);
                ug->u.a[i].a = NULL;
            }
        }
    }

    ma_ug_seq(ug, read_g, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "Writing raw unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.r_utg.noseq.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, read_g, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
    exit(0);
}


uint32_t print_untig_by_read(ma_ug_t *g, const char* name, uint32_t in, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, const char* info)
{
	uint32_t i, k, rId = (uint32_t)-1;
    if(in != (uint32_t)-1) 
    {
        rId = in;
    }
    else
    {
        for (i = 0; i < R_INF.total_reads; ++i) 
        {
            if(Get_NAME_LENGTH(R_INF, i) != strlen(name)) continue;
            if(memcmp(name, Get_NAME(R_INF, i), Get_NAME_LENGTH(R_INF, i)) == 0)
            {
                fprintf(stderr, "%s: i: %u, >%.*s\n", info, i, 
                (int)Get_NAME_LENGTH((R_INF), i), Get_NAME((R_INF), i));
                rId = i;
                break;
            }
        }
    }
    
    

    if(rId == (uint32_t)-1)
    {
        fprintf(stderr, "%s: Cannot find %s\n", info, name);
        return (uint32_t)-1;
    }

    if(sources && reverse_sources) print_specfic_read_ovlp(rId, sources, reverse_sources, info); 

    if(g != NULL)
    {
        for (i = 0; i < g->u.n; ++i) 
        {
            ma_utg_t *u = &g->u.a[i];
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                if(rId == (u->a[k]>>33))
                {
                    fprintf(stderr, "%s: %s is the %u-th read at %u-th unitig (label: %u)\n", 
                    info, name, k, i, g->g->seq[i].c);
                    return i;
                }
            }
        }
    }
	




    fprintf(stderr, "%s: %s is not at any unitig\n", info, name);
    return (uint32_t)-1;
}


void print_untig(ma_ug_t *g, uint32_t uId, const char* info, uint32_t is_print_read)
{
	uint32_t i, k;
    asg_t* nsg = g->g;

    uId = uId<<1;
    asg_arc_t *aw = asg_arc_a(nsg, uId);
    uint32_t nw = asg_arc_n(nsg, uId);

    fprintf(stderr, "\n%s: c = %u, del: %u\n", info, nsg->seq[uId>>1].c, nsg->seq[uId>>1].del);
    fprintf(stderr, "%s(%u): direction = 0...\n", info, uId>>1);
    for (i = 0; i < nw; i++)
    {
        if(aw[i].del) continue;
        fprintf(stderr, "%s: (%u) v>>1: %u, dir: %u\n", info, i, aw[i].v>>1, aw[i].v&1);
    }

    uId = uId^1;
    aw = asg_arc_a(nsg, uId);
    nw = asg_arc_n(nsg, uId);
    fprintf(stderr, "\n%s(%u): direction = 1...\n", info, uId>>1);
    for (i = 0; i < nw; i++)
    {
        if(aw[i].del) continue;
        fprintf(stderr, "%s: (%u) v>>1: %u, dir: %u\n", info, i, aw[i].v>>1, aw[i].v&1);
    }

    if(is_print_read == 0) return;


    uId = uId>>1;
    ma_utg_t *u = &g->u.a[uId];
    fprintf(stderr, "\n%s: include %u reads in total...\n", info, u->n);
    for (k = 0; k < u->n; k++)
    {
        fprintf(stderr, "%s: rId>>1: %lu, dir: %lu, name: %.*s\n", 
        info, (unsigned long)(u->a[k]>>33), (unsigned long)((u->a[k]>>32)&1),
        (int)Get_NAME_LENGTH((R_INF), (u->a[k]>>33)), Get_NAME((R_INF), (u->a[k]>>33)));
    }
}


void print_read_all(ma_ug_t *ug, const char* name, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, const char* info)
{
    fprintf(stderr, "\n\n****************************\n");
    uint32_t uId;
    uId = print_untig_by_read(ug, name, (uint32_t)-1, sources, reverse_sources, info);
    if(uId != (uint32_t)-1) print_untig(ug, uId, info, 1);
    fprintf(stderr, "****************************\n");
}


void renew_utg(ma_ug_t **ug, asg_t* read_g, kvec_asg_arc_t_warp* edge)
{
    ma_ug_t* high_level_ug = NULL;
    asg_t* nsg = (*ug)->g;
    uint32_t i;
    ma_utg_t *u;
    high_level_ug = ma_ug_gen(nsg);
    for (i = 0; i < high_level_ug->u.n; i++)
    {
        high_level_ug->g->seq[i].c = PRIMARY_LABLE;
        u = &(high_level_ug->u.a[i]);
        if(u->m == 0) continue;
        merge_unitig_content(u, (*ug), read_g, edge);
    }
    ma_ug_destroy((*ug));
    (*ug) = high_level_ug;
}



long long get_graph_statistic(asg_t *g)
{
    long long num_arc = 0;
	uint32_t n_vtx = g->n_seq * 2, v;
    
	for (v = 0; v < n_vtx; ++v) 
    {
		if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///num_arc += asg_arc_n(g, v);
        num_arc += get_real_length(g, v, NULL);
    }

    return num_arc;
}

void get_trio_labs(buf_t* b, ma_ug_t *ug, Trio_counter* flag)
{    
    flag->father_occ = flag->mother_occ = flag->ambig_occ = flag->drop_occ = 0;
    uint32_t i, v, k, rId;
    if(ug == NULL)
    {
        for (i = 0; i < b->b.n; ++i)
        {
            if(R_INF.trio_flag[b->b.a[i]>>1]==FATHER)
            {
                flag->father_occ++;
            }
            else if(R_INF.trio_flag[b->b.a[i]>>1]==MOTHER)
            {
                flag->mother_occ++;
            }
            else if(R_INF.trio_flag[b->b.a[i]>>1]==AMBIGU)
            {
                flag->ambig_occ++;
            }
            else if(R_INF.trio_flag[b->b.a[i]>>1]==DROP)
            {
                flag->drop_occ++;
            }
        }
    }
    else
    {
        ma_utg_t* u = NULL;
        for (i = 0; i < b->b.n; ++i)
        {
            v = b->b.a[i]>>1;
            u = &(ug->u.a[v]);
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                rId = u->a[k]>>33;
                if(R_INF.trio_flag[rId]==FATHER)
                {
                    flag->father_occ++;
                }
                else if(R_INF.trio_flag[rId]==MOTHER)
                {
                    flag->mother_occ++;
                }
                else if(R_INF.trio_flag[rId]==AMBIGU)
                {
                    flag->ambig_occ++;
                }
                else if(R_INF.trio_flag[rId]==DROP)
                {
                    flag->drop_occ++;
                }
            }
        }
    }

    flag->total = flag->father_occ + flag->mother_occ + flag->ambig_occ + flag->drop_occ;
}

uint32_t cal_trio_vec(buf_t* b, ma_ug_t *ug, float thres)
{
    uint32_t i, father_occ = 0, mother_occ = 0, v, k, rId;
    if(ug == NULL)
    {
        for (i = 0; i < b->b.n; ++i)
        {
            if(R_INF.trio_flag[b->b.a[i]>>1]==FATHER)
            {
                father_occ++;
            }
            else if(R_INF.trio_flag[b->b.a[i]>>1]==MOTHER)
            {
                mother_occ++;
            }
        }
    }
    else
    {
        ma_utg_t* u = NULL;
        for (i = 0; i < b->b.n; ++i)
        {
            v = b->b.a[i]>>1;
            u = &(ug->u.a[v]);
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                rId = u->a[k]>>33;
                if(R_INF.trio_flag[rId]==FATHER)
                {
                    father_occ++;
                }
                else if(R_INF.trio_flag[rId]==MOTHER)
                {
                    mother_occ++;
                }
            }
        }
    }
    
    

    if(father_occ >= thres*(father_occ+mother_occ)) return FATHER;
    if(mother_occ >= thres*(father_occ+mother_occ)) return MOTHER;
    return AMBIGU;
}

int cut_trio_tip_primary(asg_t *g, ma_ug_t *ug, uint32_t max_ext, uint32_t trio_flag, uint32_t keep_out_node,
asg_t *read_sg, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, uint32_t min_edge_length)
{
    double startTime = Get_T();
    uint32_t n_vtx = g->n_seq * 2, v, w, i, cnt = 0, tipEvaluateLen, flag, inner_flag, operation, tip_trio_flag;
    uint32_t non_trio_flag = (uint32_t)-1, nw;
    asg_arc_t *aw = NULL;
    if(trio_flag == FATHER) non_trio_flag = MOTHER;
    if(trio_flag == MOTHER) non_trio_flag = FATHER;

    buf_t b, b_0, b_1;
    
    memset(&b, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        if (g->seq[v>>1].c == CUT || g->seq[v>>1].c == CUT_DIF_HAP || g->seq[v>>1].c == TRIM) continue;
        if(get_real_length(g, (v^1), NULL) != 0) continue;
        ///cut tip of length <= max_ext
        flag = check_tip(g, v, &w, &b, max_ext);
        if(flag == LOOP || flag == LONG_TIPS) continue;
        if(keep_out_node && flag == MUL_OUTPUT) continue;
        if(ug!=NULL)
        {
            tipEvaluateLen = 0;
            for (i = 0; i < b.b.n; ++i)
            {
                tipEvaluateLen += EvaluateLen(ug->u, (b.b.a[i]>>1));
            }
            if(tipEvaluateLen > max_ext) continue;
        }

        operation = TRIM;
        ///we need to deal with MUL_INPUT and END_TIPS
        ///two operations: 1. trimming 2. cutting
        if(flag == END_TIPS)
        {
            tip_trio_flag = cal_trio_vec(&b, ug, TRIO_THRES);
            if(((tip_trio_flag == FATHER) || (tip_trio_flag == MOTHER))
                &&
                (tip_trio_flag != non_trio_flag))
            {
                operation = CUT;
            }
        }

        ///if(flag == MUL_INPUT || flag == MUL_OUTPUT)
        if(flag == MUL_INPUT)
        {
            // tip_trio_flag = cal_trio_vec(&b, ug, TRIO_THRES);
            // if(((tip_trio_flag == FATHER) || (tip_trio_flag == MOTHER))
            //     &&
            //     (tip_trio_flag != non_trio_flag))
            // {
            //     operation = CUT;
            // }
            /**
            To do lists: we can refine this function at the final step
            1. ideally, we should check if this tip comes from different haplotype with another unitig.
            That might be helpful to recover bubbles. We should keep all break points in a vector, 
            and break these types of tip. At last, recover them.
            
            2. if this tip has enough haplotype information, we should check if it does not
            come from different haplotype with another unitig. If it does not, do not trim them
            if(operation == TRIM)
            {
                get_real_length(g, b.b.a[b.b.n-1], &w);
                w = w^1;
                ;
            }
            **/
            get_real_length(g, b.b.a[b.b.n-1], &w);
            w = w^1;
            aw = asg_arc_a(g, w);
            nw = asg_arc_n(g, w);
            for (i = 0; i < nw; ++i)
            {
                if(operation == CUT) break;
                if(aw[i].del) continue;
                if(aw[i].v == (b.b.a[b.b.n-1]^1)) continue;
                inner_flag = check_different_haps(g, ug, read_sg, b.b.a[b.b.n-1]^1, aw[i].v,
                reverse_sources, &b_0, &b_1, ruIndex, min_edge_length, 1);
                if(inner_flag == NON_PLOID) operation = CUT;
            }
        }

        
        for (i = 0; i < b.b.n; ++i)
        {
            g->seq[(b.b.a[i]>>1)].c = ALTER_LABLE;
        }

        for (i = 0; i < b.b.n; ++i)
        {
            asg_seq_drop(g, (b.b.a[i]>>1));
        }
        

        if(operation == CUT)
        {
            for (i = 0; i < b.b.n; ++i)
            {
                g->seq[(b.b.a[i]>>1)].c = CUT;
            }
        }
        
		++cnt;
    }


    if (cnt > 0) asg_cleanup(g);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);


	return cnt;
}

int asg_arc_cut_trio_long_tip_primary(asg_t *g, ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t min_edge_length, float drop_ratio)
{
    // hamt note: g is auxsg
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, v_maxLen_i = (uint32_t)-1, flag, operation;
    uint32_t return_flag, n_tips;
    long long ll, v_maxLen, tmp, max_stop_nodeLen, max_stop_baseLen;

    buf_t b, b_0, b_1;
    memset(&b, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, k, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        v_maxLen = -1;
        v_maxLen_i = (uint32_t)-1;
        n_tips = 0;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            { 
                return_flag = get_unitig(g, ug, av[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, 
                &max_stop_baseLen, 1, NULL);

                /**********************for debug************************/
                // uint32_t debug_return_flag;
                // long long debug_ll;
                // debug_return_flag = get_unitig_back(g, ug, av[i].v, &convex, &debug_ll, &tmp, NULL);
                // if(debug_return_flag != return_flag || debug_ll != ll)
                // {
                //     fprintf(stderr, "ERROR\n");
                // }
                /**********************for debug************************/

                if(return_flag==LOOP) continue;
                if(return_flag==END_TIPS) n_tips++;
                
                if(v_maxLen < ll)
                {
                    v_maxLen = ll;
                    v_maxLen_i = i;
                }
            }
        }

        if(n_tips==0) continue;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                ///skip the longest way
                if(v_maxLen_i == i) continue;

                b.b.n = 0;
                if(get_unitig(g, ug, av[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, 
                &max_stop_baseLen, 1, &b)!=END_TIPS)
                {
                    continue;
                } 

                if(ll >= (v_maxLen*drop_ratio)) continue;
                n_reduced++;

                operation = TRIM;
                flag = check_different_haps(g, ug, read_sg, av[v_maxLen_i].v, av[i].v,
                reverse_sources, &b_0, &b_1, ruIndex, min_edge_length, 1);
                // #define UNAVAILABLE (uint32_t)-1
                // #define PLOID 0
                // #define NON_PLOID 1
                if(flag == NON_PLOID) operation = CUT;


                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
                }

                for (k = 0; k < b.b.n; k++)
                {
                    asg_seq_drop(g, b.b.a[k]>>1);
                }

                if(operation == CUT)
                {
                    for (k = 0; k < b.b.n; k++)
                    {
                        g->seq[b.b.a[k]>>1].c = CUT;
                    }
                }


            }
        }
    }


    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}

int asg_arc_cut_trio_long_tip_primary_complex(asg_t *g, ma_ug_t *ug, asg_t *read_sg, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t min_edge_length, float drop_ratio, uint32_t stops_threshold)
{
    double startTime = Get_T();
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, in, flag, operation;
    uint32_t return_flag, convex_i, k;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;

    buf_t b, b_0, b_1;
    memset(&b, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        if(g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if(get_real_length(g, v^1, NULL) != 1) continue;

        b.b.n = 0;   
        return_flag = get_unitig(g, ug, v^1, &convex, &ll, &tmp, &max_stop_nodeLen, 
        &max_stop_baseLen, 1, &b);

        if(return_flag != MUL_INPUT) continue;
        in = convex^1;
        get_real_length(g, convex, &convex);
        convex = convex^1;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);
        for (i = 0; i < n_convex; i++)
        {
            if(a_convex[i].del) continue;
            if(a_convex[i].v == in) break;
        }
        convex_i = i;
        ///if(convex_i == n_convex) fprintf(stderr, "ERROR\n");


        for (i = 0; i < n_convex; i++)
        {
            if(a_convex[i].del) continue;
            if(i == convex_i) continue;

            return_flag = get_unitig(g, ug, a_convex[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, 
            &max_stop_baseLen, stops_threshold, NULL);

            /**********************for debug************************/
            // uint32_t debug_return_flag, debug_convex;
            // long long debug_ll, debug_max_stop;

            // debug_return_flag = untig_detect_single_path_with_dels_n_stops(g, ug, 
            // a_convex[i].v, &debug_convex, &debug_ll, &debug_max_stop, NULL, stops_threshold);

            // if(debug_convex != convex || debug_ll != ll || debug_max_stop != max_stop_nodeLen)
            // {
            //     fprintf(stderr, "ERROR\n");
            // }
            // if(debug_return_flag != return_flag)
            // {
            //     if(
            //         !((debug_return_flag == TWO_OUTPUT && return_flag == MUL_OUTPUT)
            //         ||
            //         (debug_return_flag == TWO_INPUT && return_flag == MUL_INPUT))
            //       )
            //       {
            //           fprintf(stderr, "ERROR\n");
            //       }
            // }
            /**********************for debug************************/


            
            if(convexLen < ll*drop_ratio && max_stop_nodeLen >= ll*MAX_STOP_RATE)
            {
                n_reduced++;
                operation = TRIM;
                flag = check_different_haps(g, ug, read_sg, a_convex[convex_i].v, a_convex[i].v,
                reverse_sources, &b_0, &b_1, ruIndex, min_edge_length, stops_threshold);
                // #define UNAVAILABLE (uint32_t)-1
                // #define PLOID 0
                // #define NON_PLOID 1
                if(flag == NON_PLOID) operation = CUT;

                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
                }

                for (k = 0; k < b.b.n; k++)
                {
                    asg_seq_drop(g, b.b.a[k]>>1);
                }

                if(operation == CUT)
                {
                    for (k = 0; k < b.b.n; k++)
                    {
                        g->seq[b.b.a[k]>>1].c = CUT;
                    }
                }

                break;
            }
        }
    }


    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);

    return n_reduced;
}

void renew_longest_tip_by_drop(asg_t *g, ma_ug_t *ug, asg_arc_t *av, uint32_t nv,
long long* base_maxLen, long long* base_maxLen_i, uint32_t stops_threshold, buf_t* b, uint32_t trio_flag)
{
    if(trio_flag != FATHER && trio_flag != MOTHER) return;
    Trio_counter max, cur;
    memset(&max, 0, sizeof(Trio_counter));
    memset(&cur, 0, sizeof(Trio_counter));
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen, max_weight = 0, cur_weight = 0;
    uint32_t convex, i, return_flag, best_tip_i, best_tip_len;
    
    b->b.n = 0;
    return_flag = get_unitig(g, ug, av[(*base_maxLen_i)].v, &convex, &tmp, &ll, 
    &max_stop_nodeLen, &max_stop_baseLen, stops_threshold, b);
    if(return_flag!=END_TIPS) return;


    get_trio_labs(b, ug, &max);
    ///means this unitig might be at current haplotype
    ///if(max.drop_occ<(max.total*TRIO_DROP_THRES)) return;


    best_tip_i = (*base_maxLen_i);
    best_tip_len = (*base_maxLen);

    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if(i==(*base_maxLen_i)) continue;

        b->b.n = 0;
        return_flag = get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, 
        &max_stop_nodeLen, &max_stop_baseLen, stops_threshold, b);

        if(return_flag!=END_TIPS) return;
        ///this tip should be long enough
        if(ll<(TRIO_DROP_LENGTH_THRES*(*base_maxLen))) continue;

        get_trio_labs(b, ug, &cur);

        if(trio_flag == FATHER)
        {
            max_weight = (long long)max.father_occ - (long long)max.drop_occ - (long long)max.mother_occ;
            cur_weight = (long long)cur.father_occ - (long long)cur.drop_occ - (long long)cur.mother_occ;
            ///this unitig is very likly at another haplotype, ignore it
            if((cur.drop_occ+cur.mother_occ)>=(cur.total*TRIO_DROP_THRES)) continue;
        }

        if(trio_flag == MOTHER)
        {
            max_weight = (long long)max.mother_occ - (long long)max.drop_occ - (long long)max.father_occ;
            cur_weight = (long long)cur.mother_occ - (long long)cur.drop_occ - (long long)cur.father_occ;
            ///this unitig is very likly at another haplotype, ignore it
            if((cur.drop_occ+cur.father_occ)>=(cur.total*TRIO_DROP_THRES)) continue;
        }


        if(cur_weight > max_weight)
        {
            max = cur;
            best_tip_i = i;
            best_tip_len = ll;
        }
        else if(cur_weight == max_weight && ll>best_tip_len)
        {
            max = cur;
            best_tip_i = i;
            best_tip_len = ll;
        }
    }

    (*base_maxLen_i) = best_tip_i;
    (*base_maxLen) = best_tip_len;
}

int asg_arc_cut_trio_long_equal_tips_assembly(asg_t *g, ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex, uint32_t trio_flag)
{
    double startTime = Get_T();
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag, is_hap, n_tips, return_flag, k;
    long long ll, base_maxLen, base_maxLen_i, max_stop_nodeLen, max_stop_baseLen, tmp;
    buf_t b, b_0, b_1;
    memset(&b, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;

        base_maxLen = -1;
        base_maxLen_i = -1;
        n_tips = 0;
        is_hap = 0;

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            return_flag = get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, 
            &max_stop_baseLen, 1, NULL);

            /**********************for debug************************/
            // uint32_t debug_return_flag;
            // long long debug_ll;
            // debug_return_flag = get_unitig_back(g, ug, av[i].v, &convex, &tmp, &debug_ll, NULL);
            // if(debug_return_flag != return_flag || debug_ll != ll)
            // {
            //     fprintf(stderr, "* ERROR\n");
            // }
            /**********************for debug************************/

            if(return_flag==LOOP) continue;
            if(return_flag==END_TIPS) n_tips++;

            if(base_maxLen < ll)
            {
                base_maxLen = ll;
                base_maxLen_i = i;
            }
        }

        if(n_tips==0) continue;
        ///all the unitigs here are tips
        if(n_arc == n_tips)
        {
            renew_longest_tip_by_drop(g, ug, av, nv, &base_maxLen, 
            &base_maxLen_i, 1, &b, trio_flag);
        }

        for (i = 0; i < nv; i++)
        {
            if(i==base_maxLen_i) continue;
            if(av[i].del) continue;

            b.b.n = 0;
            return_flag = get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, 
            &max_stop_baseLen, 1, &b);

            if(return_flag != END_TIPS) continue;

            flag = check_different_haps(g, ug, read_sg, av[base_maxLen_i].v, av[i].v,
            reverse_sources, &b_0, &b_1, ruIndex, miniedgeLen, 1);

            // #define UNAVAILABLE (uint32_t)-1
            // #define PLOID 0
            // #define NON_PLOID 1
            if(flag != PLOID) continue;
            n_reduced++;

            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
            }

            for (k = 0; k < b.b.n; k++)
            {
                asg_seq_drop(g, b.b.a[k]>>1);
            }

            is_hap++;
        }

        if(is_hap > 0)
        {
            i = base_maxLen_i;
            b.b.n = 0;
            get_unitig(g, ug, av[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, 
            &max_stop_baseLen, 1, &b);

            for (k = 0; k < b.b.n; k++)
            {
                g->seq[b.b.a[k]>>1].c = HAP_LABLE;
            }
        }
    }

    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);


    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}

uint8_t if_primary_unitig(ma_utg_t* u, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex, uint8_t* r_flag)
{
    if(asm_opt.recover_atg_cov_min < 0 || asm_opt.recover_atg_cov_max < 0)
    {
        return 0;
    }
    long long R_bases = 0, C_bases = 0, C_bases_primary = 0, C_bases_alter = 0;
    long long total_C_bases = 0, total_C_bases_primary = 0;
    uint32_t available_reads = 0, k, j, rId, tn, is_Unitig;
    ma_hit_t *h;
    if(u->m == 0) return 0;
    total_C_bases = total_C_bases_primary = available_reads = 0;
    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        r_flag[rId] = 1;
    }

    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        C_bases = C_bases_primary = C_bases_alter = 0;
        R_bases = coverage_cut[rId].e - coverage_cut[rId].s;
        for (j = 0; j < (uint64_t)(sources[rId].length); j++)
        {
            h = &(sources[rId].buffer[j]);
            if(h->el != 1) continue;
            tn = Get_tn((*h));
            if(read_g->seq[tn].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }
            if(read_g->seq[tn].del == 1) continue;
            if(r_flag[tn])
            {
                C_bases_primary += Get_qe((*h)) - Get_qs((*h));
            }
            else
            {
                C_bases_alter += Get_qe((*h)) - Get_qs((*h));
            }
        }

        C_bases = C_bases_primary + C_bases_alter;
        total_C_bases += C_bases; 
        ///if(C_bases_primary < C_bases * ALTER_COV_THRES) continue;

        C_bases = C_bases/R_bases;
        if(C_bases >= asm_opt.recover_atg_cov_min && C_bases <= asm_opt.recover_atg_cov_max)
        {
            if(C_bases_primary >= (C_bases_primary + C_bases_alter) * ALTER_COV_THRES) available_reads++;
            total_C_bases_primary += C_bases_primary;
        }
    }

    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        r_flag[rId] = 0;
    }

    ///fprintf(stderr, "available_reads: %u, u->n: %u\n", available_reads, u->n);

    //if(available_reads < (u->n * 0.8) || available_reads == 0)
    if(((available_reads < (u->n * 0.8)) && (total_C_bases_primary < (total_C_bases * 0.8))) 
       || available_reads == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}



int magic_trio_phasing(asg_t *g, ma_ug_t *ug, asg_t *read_sg, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long miniedgeLen, 
R_to_U* ruIndex, uint32_t positive_flag, float drop_rate)
{
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    ma_utg_t* nsu = NULL;
    uint32_t flag = (uint32_t)-1, flag_occ, non_flag_occ, ambigious, del_node, keep_node;
    if(positive_flag == FATHER) flag = MOTHER;
    if(positive_flag == MOTHER) flag = FATHER;
    uint8_t* primary_flag = (uint8_t*)calloc(read_sg->n_seq, sizeof(uint8_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv < 2 || g->seq[v>>1].del /**|| g->seq[v>>1].c == ALTER_LABLE**/) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;
        del_node = keep_node = 0;
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nsu = &(ug->u.a[av[i].v>>1]);
            get_unitig_trio_flag(nsu, flag, &flag_occ, &non_flag_occ, &ambigious);
            ///we may need it or not
            if((flag_occ <= ((non_flag_occ+flag_occ)*drop_rate))||((flag_occ+non_flag_occ) == 0))
            {
                keep_node++;
                continue;
            }

            if(nsu->n >= 100)
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES)
                {
                    keep_node++;
                    continue;
                }
            }
            else if(nsu->n >= 50)
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES*0.5)
                {
                    keep_node++;
                    continue;
                }
            }
            else
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES*0.25)
                {
                    keep_node++;
                    continue;
                }
            }
            
            

            del_node++; 
        }

        if(keep_node == 0 || del_node == 0) continue;

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            nsu = &(ug->u.a[av[i].v>>1]);

            get_unitig_trio_flag(nsu, flag, &flag_occ, &non_flag_occ, &ambigious);
            ///we may need it or not
            if((flag_occ <= ((non_flag_occ+flag_occ)*drop_rate))||((flag_occ+non_flag_occ) == 0))
            {
                continue;
            }

            if(nsu->n >= 100)
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES)
                {
                    continue;
                }
            }
            else if(nsu->n >= 50)
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES*0.5)
                {
                    continue;
                }
            }
            else
            {
                if(flag_occ < nsu->n*DOUBLE_CHECK_THRES*0.25)
                {
                    continue;
                }
            }

            if(if_primary_unitig(nsu, read_sg, coverage_cut, sources, ruIndex, primary_flag))
            {
                continue;
            }

            g->seq[av[i].v>>1].c = ALTER_LABLE;
            asg_seq_drop(g, av[i].v>>1);
            n_reduced++;
        }
    }




    asg_cleanup(g);
    free(primary_flag);
    return n_reduced;
}

int asg_arc_cut_trio_long_equal_tips_assembly_complex(asg_t *g, ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, R_to_U* ruIndex, uint32_t stops_threshold)
{
    double startTime = Get_T();
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0, convex, in, flag;
    uint32_t return_flag, convex_i, k;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;

    buf_t b, b_0, b_1;
    memset(&b, 0, sizeof(buf_t));
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i;
        if (g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
        ///tip
        if (get_real_length(g, v, NULL) != 0) continue;
        if (get_real_length(g, v^1, NULL) != 1) continue;

        b.b.n = 0;   
        return_flag = get_unitig(g, ug, v^1, &convex, &tmp, &ll, &max_stop_nodeLen, 
        &max_stop_baseLen, 1, &b);

        if(return_flag != MUL_INPUT) continue;
        in = convex^1;
        get_real_length(g, convex, &convex);
        convex = convex^1;
        uint32_t n_convex = asg_arc_n(g, convex), convexLen = ll;
        asg_arc_t *a_convex = asg_arc_a(g, convex);
        for (i = 0; i < n_convex; i++)
        {
            if(a_convex[i].del) continue;
            if(a_convex[i].v == in) break;
        }
        convex_i = i;
        ///if(convex_i == n_convex) fprintf(stderr, "ERROR1\n");

        for (i = 0; i < n_convex; i++)
        {
            if(a_convex[i].del) continue;
            if(i == convex_i) continue;

            return_flag = get_unitig(g, ug, a_convex[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, 
            &max_stop_baseLen, stops_threshold, NULL);


            /**********************for debug************************/
            // uint32_t debug_return_flag, debug_convex;
            // long long debug_ll, debug_max_stop;

            // debug_return_flag = untig_detect_single_path_with_dels_contigLen_complex(g,  
            // a_convex[i].v, &debug_convex, &debug_ll, &debug_max_stop, NULL, stops_threshold);

            // if(debug_convex != convex || debug_ll != ll || debug_max_stop != max_stop_baseLen)
            // {
            //     fprintf(stderr, "ERROR2\n");
            // }
            // if(debug_return_flag != return_flag)
            // {
            //     if(
            //         !((debug_return_flag == TWO_OUTPUT && return_flag == MUL_OUTPUT)
            //         ||
            //         (debug_return_flag == TWO_INPUT && return_flag == MUL_INPUT))
            //       )
            //       {
            //           fprintf(stderr, "ERROR3\n");
            //       }
            // }
            /**********************for debug************************/

            if(ll>convexLen  && max_stop_baseLen>=ll*MAX_STOP_RATE)
            {
                flag = check_different_haps(g, ug, read_sg, a_convex[convex_i].v, a_convex[i].v,
                reverse_sources, &b_0, &b_1, ruIndex, miniedgeLen, stops_threshold);

                // #define UNAVAILABLE (uint32_t)-1
                // #define PLOID 0
                // #define NON_PLOID 1
                if(flag != PLOID) continue;
                n_reduced++;

                for (k = 0; k < b.b.n; k++)
                {
                    g->seq[b.b.a[k]>>1].c = ALTER_LABLE;
                }

                for (k = 0; k < b.b.n; k++)
                {
                    asg_seq_drop(g, b.b.a[k]>>1);
                }


                ///lable the primary one
                b_0.b.n = 0; 
                get_unitig(g, ug, a_convex[i].v, &convex, &tmp, &ll, &max_stop_nodeLen, 
                &max_stop_baseLen, stops_threshold, &b_0);
                for (k = 0; k < b_0.b.n; k++)
                {
                    g->seq[b_0.b.a[k]>>1].c = HAP_LABLE;
                }

                break;
            }
        }

    }

    asg_cleanup(g);
    asg_symm(g);
    free(b.b.a);
    free(b_0.b.a);
    free(b_1.b.a);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}


int detect_chimeric_by_topo(asg_t *g, ma_ug_t *ug, asg_t *read_sg, 
ma_hit_t_alloc* reverse_sources, long long miniedgeLen, uint32_t stops_threshold, float drop_rate,
R_to_U* ruIndex)
{
    double startTime = Get_T();
    uint32_t i, k, v_i, v_beg, v_end, selfLen, w1, w2, wv, nw, n_vtx = g->n_seq * 2, n_reduced = 0, convex, convex_T, read_num;
    asg_arc_t *aw;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;

    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    for (v_i = 0; v_i < n_vtx; ++v_i) 
    {
        v_beg = v_i;
        ///some node could be deleted
        if (g->seq[v_beg>>1].del || g->seq[v_beg>>1].c == ALTER_LABLE) continue;
        if(get_real_length(g, v_beg, NULL) != 1) continue;

        get_real_length(g, v_beg, &w1);
        if(get_real_length(g, w1^1, NULL)<=1) continue;

        if(get_unitig(g, ug, v_beg^1, &v_end, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
        1, NULL)==LOOP)
        {
            continue;
        }
        selfLen = ll;
        if(get_real_length(g, v_end, NULL) != 1) continue;


        get_real_length(g, v_end, &w2);
        if(get_real_length(g, w2^1, NULL)<=1) continue;

        get_unitig(g, ug, w1, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
        stops_threshold, NULL);
        if(ll*drop_rate < selfLen) continue;
        if(ll*0.4 > max_stop_nodeLen) continue;

        get_unitig(g, ug, w2, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
        stops_threshold, NULL);
        if(ll*drop_rate < selfLen) continue;
        if(ll*0.4 > max_stop_nodeLen) continue;

        stops_threshold++;


        w1 = w1^1;wv=v_beg^1;aw = asg_arc_a(g, w1); nw = asg_arc_n(g, w1);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            get_unitig(g, ug, aw[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
            stops_threshold, NULL);
            if(ll*drop_rate < selfLen) break;
            if(ll*0.4 > max_stop_nodeLen) continue;

            if((aw[i].v>>1)==(v_beg>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v_beg>>1)) continue;

            b_0.b.n = 0;
            get_unitig(g, ug, aw[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
            stops_threshold, &b_0);
            ///if(convex == convex_T) break;
            for (k = 0; k < b_0.b.n; k++)
            {
                if((b_0.b.a[k]>>1) == (convex_T>>1)) break;
            }
            if(k != b_0.b.n) break;
            
            if(check_different_haps(g, ug, read_sg, wv, aw[i].v, reverse_sources, &b_0, &b_1, 
            ruIndex, miniedgeLen, stops_threshold)==PLOID)
            {
                break;
            }
        }
        if(i!=nw) continue;





        w2 = w2^1;wv=v_end^1;aw = asg_arc_a(g, w2); nw = asg_arc_n(g, w2);convex_T = (uint32_t)-1;
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;

            
            get_unitig(g, ug, aw[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
            stops_threshold, NULL);
            if(ll*drop_rate < selfLen) break;
            if(ll*0.4 > max_stop_nodeLen) continue;

            if((aw[i].v>>1) == (v_end>>1))
            {
                if(convex_T != (uint32_t)-1) break;
                convex_T = convex;
            }
        }
        if(i!=nw) continue;

        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if((aw[i].v>>1) == (v_end>>1)) continue;

            b_0.b.n = 0;
            get_unitig(g, ug, aw[i].v, &convex, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
            stops_threshold, &b_0);
            ///if(convex == convex_T) break;
            for (k = 0; k < b_0.b.n; k++)
            {
                if((b_0.b.a[k]>>1) == (convex_T>>1)) break;
            }
            if(k != b_0.b.n) break;

            if(check_different_haps(g, ug, read_sg, wv, aw[i].v, reverse_sources, &b_0, &b_1, 
            ruIndex, miniedgeLen, stops_threshold)==PLOID)
            {
                break;
            }
        }
        if(i!=nw) continue;

        n_reduced++;
        b_0.b.n = 0;
        read_num = 0;
        ///get_long_tip_length(g, &(ug->u), v_beg^1, &v_end, &b_0);
        get_unitig(g, ug, v_beg^1, &v_end, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 1, &b_0);
        for (k = 0; k < b_0.b.n; k++)
        {
            g->seq[b_0.b.a[k]>>1].c = ALTER_LABLE;
            read_num += ug->u.a[b_0.b.a[k]>>1].n;
        }

        for (k = 0; k < b_0.b.n; k++)
        {
            asg_seq_drop(g, b_0.b.a[k]>>1);
        }

        if(read_num <= CHIMERIC_TRIM_THRES)
        {
            for (k = 0; k < b_0.b.n; k++)
            {
                g->seq[b_0.b.a[k]>>1].c = CUT;
            }
        }

    }

    asg_cleanup(g);
    free(b_0.b.a);
    free(b_1.b.a);
    
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d long tips\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    

    return n_reduced;
}

uint32_t cmp_untig_graph(ma_ug_t *src, ma_ug_t *dest)
{
    uint32_t nv, nw;
	asg_arc_t *av = NULL, *aw = NULL;
    uint32_t v, i, k, n_vtx = dest->g->n_seq*2;
    if(src->g->n_seq!=dest->g->n_seq)
    {
        fprintf(stderr, "src->g->n_seq: %u, dest->g->n_seq: %u\n", src->g->n_seq, dest->g->n_seq);
        ///return 1;
    } 
    if(src->g->r_seq!=dest->g->r_seq)
    {
        fprintf(stderr, "src->g->r_seq: %u, dest->g->r_seq: %u\n", src->g->r_seq, dest->g->r_seq);
        ///return 1;
    } 
    if(src->g->is_srt!=dest->g->is_srt)
    {
        fprintf(stderr, "src->g->is_srt: %u, dest->g->is_srt: %u\n", src->g->is_srt, dest->g->is_srt);
        ///return 1;
    } 
    if(src->g->is_symm!=dest->g->is_symm)
    {
        fprintf(stderr, "src->g->is_symm: %u, dest->g->is_symm: %u\n", src->g->is_symm, dest->g->is_symm);
        ///return 1;
    } 
    ///if(memcmp(src->g->seq, dest->g->seq, sizeof(asg_seq_t)*src->g->n_seq)!=0) return 1;

    n_vtx = dest->g->n_seq;
    for (v = 0; v < n_vtx; v++)
    {
        if(src->g->seq[v].c != dest->g->seq[v].c)
        {
            fprintf(stderr, "src->g->seq[%u].c: %u, dest->g->seq[%u].c: %u\n", 
            v, src->g->seq[v].c, v, dest->g->seq[v].c);
        }

        if(src->g->seq[v].del != dest->g->seq[v].del)
        {
            fprintf(stderr, "src->g->seq[%u].del: %u, dest->g->seq[%u].del: %u\n", 
            v, src->g->seq[v].del, v, dest->g->seq[v].del);
        }

        if(src->g->seq[v].len != dest->g->seq[v].len)
        {
            fprintf(stderr, "src->g->seq[%u].len: %u, dest->g->seq[%u].len: %u\n", 
            v, src->g->seq[v].len, v, dest->g->seq[v].len);
        }
    }
    

    n_vtx = dest->g->n_seq*2;
    for (v = 0; v < n_vtx; v++)
    {
        nv = asg_arc_n(src->g, v);
	    av = asg_arc_a(src->g, v);
        nw = asg_arc_n(dest->g, v);
        aw = asg_arc_a(dest->g, v);
        if(get_real_length(src->g, v, NULL) != get_real_length(dest->g, v, NULL))
        {
            fprintf(stderr, "v>>1: %u, v&1: %u, get_real_length(src->g, v, NULL): %u, get_real_length(dest->g, v, NULL): %u\n", 
            v>>1, v&1, get_real_length(src->g, v, NULL), get_real_length(dest->g, v, NULL));
            ///return 1;
        } 

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(av[i].v == aw[k].v) break;
            }

            if(k == nw) return 1;
        }
    }

    n_vtx = dest->g->n_seq;
    ma_utg_t *src_u = NULL, *dest_u = NULL;
    for (v = 0; v < n_vtx; v++)
    {
        src_u = &(src->u.a[v]);
        dest_u = &(dest->u.a[v]);
        if(src_u->circ!=dest_u->circ) return 1;
        if(src_u->end!=dest_u->end) return 1;
        if(src_u->len!=dest_u->len) return 1;
        if(src_u->n!=dest_u->n) return 1;
        if(src_u->start!=dest_u->start) return 1;
        if(memcmp(src_u->a, dest_u->a, 8*src_u->n)!=0) return 1;
    }

    return 0;
}

ma_ug_t* copy_untig_graph(ma_ug_t *src)
{
    ma_ug_t *ug = NULL;
    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = asg_init();
    uint32_t v;

    ug->g->n_F_seq = src->g->n_F_seq;
    ug->g->r_seq = src->g->r_seq;
    ug->g->m_seq = ug->g->n_seq = src->g->n_seq;
    ug->g->seq = (asg_seq_t*)malloc(ug->g->n_seq * sizeof(asg_seq_t));
    memcpy(ug->g->seq, src->g->seq, sizeof(asg_seq_t)*ug->g->n_seq);

    ug->g->m_arc = ug->g->n_arc = src->g->n_arc;
    ug->g->arc = (asg_arc_t*)malloc(ug->g->n_arc*sizeof(asg_arc_t));
    memcpy(ug->g->arc, src->g->arc, sizeof(asg_arc_t)*ug->g->n_arc);

    ug->g->is_srt = src->g->is_srt;
    ug->g->is_symm = src->g->is_symm;
    ug->g->idx = (uint64_t*)malloc(ug->g->n_seq*2*8);
    memcpy(ug->g->idx, src->g->idx, ug->g->n_seq*2*8); 
    asg_cleanup(ug->g);



    ug->u.m = ug->u.n = src->u.n;
    ug->u.a = (ma_utg_t*)malloc(sizeof(ma_utg_t)*ug->u.n);

    ma_utg_t *src_u = NULL, *ug_u = NULL;
    for (v = 0; v < ug->u.n; v++)
    {
        src_u = &(src->u.a[v]);
        ug_u = &(ug->u.a[v]);
        (*ug_u) = (*src_u);
        ug_u->s = NULL;
        ug_u->a = NULL;
        ug_u->m = ug_u->n = src_u->n;

        ug_u->a = (uint64_t*)malloc(8 * ug_u->n);
        memcpy(ug_u->a, src_u->a, 8*ug_u->n);
    }

    /**
    if(cmp_untig_graph(src, ug)) fprintf(stderr, "ERROR\n");
    ma_ug_destroy(ug);
    return NULL;
    **/
    return ug;
}

void clean_trio_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long bubble_dist, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean, int just_bubble_pop, 
float drop_ratio, uint32_t trio_flag, float trio_drop_rate)
{
    asg_t *g = ug->g;
    uint32_t is_first = 1;

    redo:
    ///print_untig((ug), 61955, "i-0:", 0);

    asg_pop_bubble_primary_trio(ug, bubble_dist, trio_flag, DROP);
    untig_asg_arc_simple_large_bubbles_trio(ug, read_g, reverse_sources, 2, ruIndex, trio_flag, DROP); 
    magic_trio_phasing(g, ug, read_g, coverage_cut, sources, reverse_sources, 2, ruIndex, trio_flag, trio_drop_rate);    
    ///drop_semi_circle(ug, g, read_g, reverse_sources, ruIndex);
    /**********debug**********/
    if(just_bubble_pop == 0)
    {
        ///asg_cut_tip_primary(g, ug, tipsLen);
        cut_trio_tip_primary(g, ug, tipsLen, trio_flag, 0, read_g, reverse_sources, ruIndex, 
        2);
    }
    /**********debug**********/
    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        ///need consider tangles
        ///asg_pop_bubble_primary(g, bubble_dist);
        asg_pop_bubble_primary_trio(ug, bubble_dist, trio_flag, DROP);
        /**********debug**********/
        if(just_bubble_pop == 0)
        {
            ///need consider tangles
            asg_arc_cut_trio_long_tip_primary(g, ug, read_g, reverse_sources, ruIndex, 2, tip_drop_ratio);            
            asg_arc_cut_trio_long_equal_tips_assembly(g, ug, read_g, reverse_sources, 2, ruIndex, trio_flag);            
            asg_arc_cut_trio_long_tip_primary_complex(g, ug, read_g, reverse_sources, ruIndex,
            2, tip_drop_ratio, stops_threshold);
            asg_arc_cut_trio_long_equal_tips_assembly_complex(g, ug, read_g, reverse_sources, 
            2, ruIndex, stops_threshold);
            ///print_debug_gfa(read_g, ug, coverage_cut, "debug_chimeric", sources, ruIndex, asm_opt.max_hang_Len, asm_opt.min_overlap_Len);
            detect_chimeric_by_topo(g, ug, read_g, reverse_sources, 2, stops_threshold, chimeric_rate,
            ruIndex);
            ///need consider tangles
            ///note we need both the read graph and the untig graph
        }
        /**********debug**********/
        cur_cons = get_graph_statistic(g);
    }
    untig_asg_arc_simple_large_bubbles_trio(ug, read_g, reverse_sources, 2, ruIndex, trio_flag, DROP);

    if(just_bubble_pop == 0)
    {   
        ///asg_cut_tip_primary(g, ug, tipsLen);
        cut_trio_tip_primary(g, ug, tipsLen, trio_flag, 0, read_g, reverse_sources, ruIndex, 
        2);
    }

    resolve_tangles(ug, read_g, reverse_sources, 20, 100, 0.05, 0.2, ruIndex, trio_flag, drop_ratio);    
    drop_semi_circle(ug, g, read_g, reverse_sources, ruIndex);
    all_to_all_deduplicate(ug, read_g, coverage_cut, sources, trio_flag, trio_drop_rate, reverse_sources, ruIndex, DOUBLE_CHECK_THRES);

    if(is_first)
    {
        is_first = 0;
        unitig_arc_del_short_diploid_by_length(ug->g, drop_ratio);
        goto redo;
    }    
}


void clean_primary_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean, int just_bubble_pop, 
float drop_ratio)
{
    #define T_ROUND 2
    asg_t *g = ug->g;
    int round = T_ROUND;
    

    redo:
    asg_pop_bubble_primary_trio(ug, bubble_dist, (uint32_t)-1, DROP);
    untig_asg_arc_simple_large_bubbles_trio(ug, read_g, reverse_sources, 2, ruIndex, (uint32_t)-1, DROP); 

    if(just_bubble_pop == 0)
    {
        cut_trio_tip_primary(g, ug, tipsLen, (uint32_t)-1, 0, read_g, reverse_sources, ruIndex, 
        2);
    }

    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        asg_pop_bubble_primary_trio(ug, bubble_dist, (uint32_t)-1, DROP);

        if(just_bubble_pop == 0)
        {
            ///need consider tangles
            // (trim/cut tips except the longest one at each spot)
            asg_arc_cut_trio_long_tip_primary(g, ug, read_g, reverse_sources, ruIndex, 
            2, tip_drop_ratio);

            asg_arc_cut_trio_long_equal_tips_assembly(g, ug, read_g, reverse_sources, 2, ruIndex, (uint32_t)-1);
            asg_arc_cut_trio_long_tip_primary_complex(g, ug, read_g, reverse_sources, ruIndex,
            2, tip_drop_ratio, stops_threshold);
            asg_arc_cut_trio_long_equal_tips_assembly_complex(g, ug, read_g, reverse_sources, 
            2, ruIndex, stops_threshold);
            detect_chimeric_by_topo(g, ug, read_g, reverse_sources, 2, stops_threshold, chimeric_rate,
            ruIndex);
            
            if(round != T_ROUND)
            {
                unitig_arc_del_short_diploid_by_length_topo(g, ug, drop_ratio, asm_opt.max_short_tip, 
                reverse_sources, 0, 1);
            }
        }
        cur_cons = get_graph_statistic(g);
    }
    untig_asg_arc_simple_large_bubbles_trio(ug, read_g, reverse_sources, 2, ruIndex, (uint32_t)-1, DROP);

    if(just_bubble_pop == 0)
    {   
        cut_trio_tip_primary(g, ug, tipsLen, (uint32_t)-1, 0, read_g, reverse_sources, ruIndex, 
        2);
    }


    
    resolve_tangles(ug, read_g, reverse_sources, 20, 100, 0.05, 0.2, ruIndex, (uint32_t)-1, drop_ratio);
    drop_semi_circle(ug, g, read_g, reverse_sources, ruIndex);

    unitig_arc_del_short_diploid_by_length_topo(g, ug, drop_ratio, asm_opt.max_short_tip, reverse_sources, 0, 1);

    if(round > 0)
    {
        if(round != T_ROUND)
        {
            unitig_arc_del_short_diploid_by_length(ug->g, drop_ratio);        
        }
        round--;
        goto redo;
    }    
}

void set_drop_trio_flag(ma_ug_t *ug)
{
    ma_utg_t* u = NULL;
    asg_t* nsg = ug->g;
    uint32_t k, rId;
    uint32_t v, n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        u = &(ug->u.a[v]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            if(R_INF.trio_flag[rId] != AMBIGU) continue;
            R_INF.trio_flag[rId] = DROP;
        }
    }
}



void update_unitig_graph(ma_ug_t* ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, 
uint8_t is_final_check, float double_check_rate, uint8_t flag, float drop_rate)
{
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq, k, rId, flag_occ, non_flag_occ, hap_label_occ, n_reduce = 1;
    ma_utg_t *u;
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));

    drop_semi_circle(ug, nsg, read_g, reverse_sources, ruIndex);

    while (n_reduce)
    {
        n_reduce = 0;
        n_vtx = nsg->n_seq;
        for (v = 0; v < n_vtx; ++v) 
        {
            if (nsg->seq[v].del) continue;
            u = &((ug)->u.a[v]);
            if(u->m == 0) continue;
            if((get_real_length(nsg, v<<1, NULL)!=0) 
                && (get_real_length(nsg, ((v<<1)^1), NULL)!=0)) continue;

            flag_occ = non_flag_occ = hap_label_occ = 0;
            for (k = 0; k < u->n; k++)
            {
                rId = u->a[k]>>33;
                if(read_g->seq[rId].c == HAP_LABLE) hap_label_occ++;
                if(R_INF.trio_flag[rId] == AMBIGU) continue;
                if(R_INF.trio_flag[rId] == DROP) continue;
                if(R_INF.trio_flag[rId] == flag) flag_occ++;
                if(R_INF.trio_flag[rId] != flag) non_flag_occ++;
            }

            if(hap_label_occ == u->n) continue;
            ///if(is_double_check && non_flag_occ < u->n*DOUBLE_CHECK_THRES) continue;
            ///if(is_double_check && non_flag_occ < u->n*double_check_rate) continue;
            if(is_final_check)
            {
                if(non_flag_occ < u->n*double_check_rate) continue;
            }
            else
            {
                if(u->n >= 100)
                {
                    if(non_flag_occ < u->n*double_check_rate)
                    {
                        continue;
                    }
                }
                else if(u->n >= 50)
                {
                    if(non_flag_occ < u->n*double_check_rate*0.5)
                    {
                        continue;
                    }
                }
                else
                {
                    if(non_flag_occ < u->n*double_check_rate*0.25)
                    {
                        continue;
                    }
                }
            }
            

            if(non_flag_occ > ((non_flag_occ+flag_occ)*drop_rate))
            {
                if(if_primary_unitig(u, read_g, coverage_cut, sources, ruIndex, primary_flag))
                {
                    continue;
                }                         
                if(u->m != 0)
                {
                    u->circ = u->end = u->len = u->m = u->n = u->start = 0;
                    free(u->a);
                    u->a = NULL;
                }
                asg_seq_del(nsg, v);
                n_reduce++;
            }
        }
    }

    drop_semi_circle(ug, nsg, read_g, reverse_sources, ruIndex);
    asg_cleanup(nsg);
    free(primary_flag);
}




void get_candidate_uids(asg_t* nsg, ma_utg_t* nsu, kvec_t_u64_warp* u_vecs, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    uint32_t k, j, rId, uId, is_Unitig, m;
    uint64_t pre;
    u_vecs->a.n = 0;
    if(nsu->m == 0) return;
    for (k = 0; k < nsu->n; k++)
    {
        rId = nsu->a[k]>>33;
        for (j = 0; j < reverse_sources[rId].length; j++)
        {
            get_R_to_U(ruIndex, Get_tn(reverse_sources[rId].buffer[j]), &uId, &is_Unitig);
            if(uId==(uint32_t)-1) continue;
            ///contained read
            if(is_Unitig == 0) get_R_to_U(ruIndex, uId, &uId, &is_Unitig);
            if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
            ///need this line
            if(nsg->seq[uId].c == ALTER_LABLE) continue;
            pre = uId; 
            pre = pre | (0x100000000);
            kv_push(uint64_t, u_vecs->a, pre);
        }
    }
    if(u_vecs->a.n == 0) return;
    radix_sort_arch64(u_vecs->a.a, u_vecs->a.a + u_vecs->a.n);
    for (m = 0, k = 1; k < u_vecs->a.n; k++)
    {
        if(u_vecs->a.a[k] == u_vecs->a.a[k-1])
        {
            u_vecs->a.a[m] += (0x100000000);
        }
        else
        {
            m++;
            u_vecs->a.a[m] = u_vecs->a.a[k];
        }
    }

    u_vecs->a.n = m + 1;
}

uint32_t unitig_simi(uint32_t x, uint32_t y, ma_ug_t* ug, ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex)
{
    uint32_t k, j, uId, tn, is_Unitig, rId, ref_unitig, min_count, max_count;
    ma_utg_t *nsu_x = NULL, *nsu_y = NULL, *nsu_query = NULL;
    nsu_x = &(ug->u.a[x]);
    nsu_y = &(ug->u.a[y]);

    if(nsu_x->n == 0 || nsu_y->n == 0) return UNAVAILABLE;

    if(nsu_x->n <= nsu_y->n)
    {
        nsu_query = nsu_x;
        ref_unitig = y;
    }
    else
    {
        nsu_query = nsu_y;
        ref_unitig = x;
    }
     
    min_count = max_count = 0;
    for (k = 0; k < nsu_query->n; k++)
    {
        rId = nsu_query->a[k]>>33;
        if(reverse_sources[rId].length >= 0) min_count++;

        for (j = 0; j < reverse_sources[rId].length; j++)
        {
            tn = Get_tn(reverse_sources[rId].buffer[j]);
            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId==(uint32_t)-1) continue;
            ///contained read
            if(is_Unitig == 0) get_R_to_U(ruIndex, uId, &uId, &is_Unitig);
            if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
            if(uId == ref_unitig)
            {
                max_count++;
                break;
            }
        }
    }
    
    if(min_count == 0) return UNAVAILABLE;
	if(max_count > min_count*DIFF_HAP_RATE) return PLOID;
	return NON_PLOID;
}

void get_unitig_trio_flag(ma_utg_t* nsu, uint32_t flag, uint32_t* require, 
uint32_t* non_require, uint32_t* ambigious)
{
    uint32_t k, rId;
    (*require) = (*non_require) = (*ambigious) = 0;
    for (k = 0; k < nsu->n; k++)
    {
        rId = nsu->a[k]>>33;
        if(R_INF.trio_flag[rId] == AMBIGU || R_INF.trio_flag[rId] == DROP)
        {
            (*ambigious)++;
            continue;
        } 
        if(R_INF.trio_flag[rId] == flag)
        {
            (*require)++;
            continue;
        } 
        if(R_INF.trio_flag[rId] != flag)
        {
            (*non_require)++;
            continue;
        } 
    }
}


///note: to use this function, don't renew unitig graph!!!!!!!!!
void all_to_all_deduplicate(ma_ug_t* ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, uint8_t postive_flag, float drop_rate, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, float double_check_rate)
{
    

    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    uint32_t n_vtx, v, k, is_Unitig, uId, rId, convex, flag_occ, non_flag_occ, ambigious, /**is_ambigious,**/ flag, n_reduce = 1;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    nsg = ug->g;
    if(postive_flag == FATHER) flag = MOTHER;
    if(postive_flag == MOTHER) flag = FATHER;
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
     
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c==ALTER_LABLE) continue;

        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }

    n_reduce = 1;
    while (n_reduce)
    {
        n_reduce = 0;
        n_vtx = nsg->n_seq*2;
        for (v = 0; v < n_vtx; ++v) 
        {
            
            uId = v>>1;
            nsu = &(ug->u.a[uId]);
            if(nsu->m == 0) continue;
            if(nsg->seq[uId].del) continue;
            if(nsg->seq[uId].c == ALTER_LABLE) continue;

            if(get_real_length(nsg, v^1, NULL) == 1)
            {
                get_real_length(nsg, v^1, &convex);
                if(get_real_length(nsg, convex^1, NULL) == 1) continue;
            }

            get_unitig_trio_flag(nsu, flag, &flag_occ, &non_flag_occ, &ambigious);
            /**
            is_ambigious = 0;
            if((flag_occ+non_flag_occ)==0 && ambigious > 0)
            {
                is_ambigious = 1;
            }
            else**/
            {
                ///we may need it or not
                if(flag_occ <= ((non_flag_occ+flag_occ)*drop_rate)) continue;
                if((flag_occ+non_flag_occ) == 0) continue;
                if(nsu->n >= 100)
                {
                    if(flag_occ < nsu->n*double_check_rate)
                    {
                        continue;
                    }
                }
                else if(nsu->n >= 50)
                {
                    if(flag_occ < nsu->n*double_check_rate*0.5)
                    {
                        continue;
                    }
                }
                else
                {
                    if(flag_occ < nsu->n*double_check_rate*0.25)
                    {
                        continue;
                    }
                }
            }
            

            get_candidate_uids(nsg, nsu, &u_vecs, reverse_sources, ruIndex);


            if(u_vecs.a.n == 0) continue;

            for (k = 0; k < u_vecs.a.n; k++)
            {
                /**
                if(is_ambigious)
                {
                    get_unitig_trio_flag(&(ug->u.a[(uint32_t)(u_vecs.a.a[k])]), postive_flag, 
                    &flag_occ, &non_flag_occ, &ambigious);
                    if(flag_occ <= ((non_flag_occ+flag_occ)*drop_rate)) continue;
                    if((flag_occ+non_flag_occ) == 0) continue;
                }**/
                if(unitig_simi(uId, (uint32_t)(u_vecs.a.a[k]), ug, reverse_sources, ruIndex)==PLOID)
                {
                    break;
                }
            }

            if(k != u_vecs.a.n)
            {
                if(if_primary_unitig(nsu, read_g, coverage_cut, sources, ruIndex, primary_flag))
                {
                    continue;
                }
                if(nsu->m != 0)
                {
                    nsu->circ = nsu->end = nsu->len = nsu->m = nsu->n = nsu->start = 0;
                    free(nsu->a);
                    nsu->a = NULL;
                }
                asg_seq_del(nsg, uId);
                n_reduce++;
            }
        }
    }

    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    kv_destroy(u_vecs.a);
    free(primary_flag);
}

void delete_useless_nodes(ma_ug_t **ug)
{
    asg_t* nsg = (*ug)->g;
    uint32_t v, n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE)
        {
            asg_seq_del(nsg, v);
            if((*ug)->u.a[v].m!=0)
            {
                (*ug)->u.a[v].m = (*ug)->u.a[v].n = 0;
                free((*ug)->u.a[v].a);
                (*ug)->u.a[v].a = NULL;
            }

            continue;
        }

        //note: after cleaning, some cirle might be gone, or we have some new circles
        //so need to renew .circ
        /**
        if(get_unitig(nsg, NULL, (v<<1), &convex, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, NULL)==LOOP)
        {
            (*ug)->u.a[v].circ = 1;
            (*ug)->u.a[v].start = (*ug)->u.a[v].end = UINT32_MAX;
        }
        else
        {
            (*ug)->u.a[v].circ = 0;

            (*ug)->u.a[v].start = (*ug)->u.a[v].a[0]>>32;
            (*ug)->u.a[v].end = ((*ug)->u.a[v].a[(*ug)->u.a[v].n-1]>>32)^1;
        }
        **/        
    }

    asg_cleanup(nsg);
}




void delete_useless_trio_nodes(ma_ug_t **ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex)
{
    asg_t* nsg = (*ug)->g;
    uint32_t v, n_vtx = nsg->n_seq;
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE &&
          (if_primary_unitig(&((*ug)->u.a[v]), read_g, coverage_cut, sources, 
          ruIndex, primary_flag) == 0))
        {
            asg_seq_del(nsg, v);
            if((*ug)->u.a[v].m!=0)
            {
                (*ug)->u.a[v].m = (*ug)->u.a[v].n = 0;
                free((*ug)->u.a[v].a);
                (*ug)->u.a[v].a = NULL;
            }

            continue;
        }

        //note: after cleaning, some cirle might be gone, or we have some new circles
        //so need to renew .circ
        /**
        if(get_unitig(nsg, NULL, (v<<1), &convex, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, NULL)==LOOP)
        {
            (*ug)->u.a[v].circ = 1;
            (*ug)->u.a[v].start = (*ug)->u.a[v].end = UINT32_MAX;
        }
        else
        {
            (*ug)->u.a[v].circ = 0;

            (*ug)->u.a[v].start = (*ug)->u.a[v].a[0]>>32;
            (*ug)->u.a[v].end = ((*ug)->u.a[v].a[(*ug)->u.a[v].n-1]>>32)^1;
        }
        **/        
    }

    asg_cleanup(nsg);
    free(primary_flag);
}





void drop_semi_circle(ma_ug_t *ug, asg_t* nsg, asg_t* read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    // note: nst is ug->g aka the auxsg
    // FUNC
    //    break soap bubble
    // NOTE
    //    break the structure at an arbitrary point, not necessarily at the "stem"
    uint32_t v, n_vtx = nsg->n_seq*2, convex_f, convex_b, i, nv;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    asg_arc_t *av = NULL;
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    for (v = 0; v < n_vtx; ++v) 
    {
        // require v to not be a tip
        if(get_real_length(nsg, v, NULL) == 0) continue;
        if(get_real_length(nsg, v^1, NULL) == 0) continue;

        // check left side
        if(get_unitig(nsg, ug, v^1, &convex_b, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
                    1, NULL)!=MUL_INPUT)
        {
            continue;
        }
        get_real_length(nsg, convex_b, &convex_b);  // convext_b now is the branching node on the left

        // check right side
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;

            if(get_unitig(nsg, ug, av[i].v, &convex_f, &ll, &tmp, &max_stop_nodeLen, &max_stop_baseLen, 
                    1, NULL)!=MUL_INPUT)
            {
                continue;
            }
            get_real_length(nsg, convex_f, &convex_f);  // convex_f is the branching node on the right
            if(convex_f != convex_b) continue;
            // if we reach here, v is in a "loop" (not a circle) whose "stem" or handle is convex_f (aka convex_b)
            if(check_different_haps(nsg, ug, read_g, v^1, av[i].v, 
                                    reverse_sources, &b_0, &b_1, ruIndex, 2, 1) == PLOID)
            {
                av[i].del = 1;
		        asg_arc_del(nsg, av[i].v^1, v^1, 1);
            }
        }
        
    }

    free(b_0.b.a);
    free(b_1.b.a);

}

void update_hap_label(ma_ug_t *ug, asg_t* read_g)
{
    uint32_t v, n_vtx, k;
    uint64_t rId;

    if(ug == NULL)
    {
        n_vtx = read_g->n_seq;
        for (v = 0; v < n_vtx; ++v) 
        {
            if(read_g->seq[v].del) continue;
            if(read_g->seq[v].c != HAP_LABLE) continue;
            read_g->seq[v].c = PRIMARY_LABLE;
        }
        return;
    }


    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    ma_utg_t* u = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c != HAP_LABLE) continue;
        u = &((ug)->u.a[v]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            read_g->seq[rId].c = HAP_LABLE;
        }
    }
}


uint8_t* get_utg_attributes(ma_ug_t *ug, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex)
{
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq, k, j, rId, available_reads = 0, tn, is_Unitig;
    ma_utg_t* u = NULL;
    ma_hit_t *h;
    long long R_bases = 0, C_bases = 0, C_bases_primary = 0, C_bases_alter = 0;
    uint8_t* r_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
    uint8_t* u_flag = (uint8_t*)calloc(n_vtx, sizeof(uint8_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        u = &(ug->u.a[v]);
        if(u->m == 0) continue;
        available_reads = 0;
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            r_flag[rId] = 1;
        }

        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            C_bases = C_bases_primary = C_bases_alter = 0;
            R_bases = coverage_cut[rId].e - coverage_cut[rId].s;
            for (j = 0; j < (uint64_t)(sources[rId].length); j++)
            {
                h = &(sources[rId].buffer[j]);
                if(h->el != 1) continue;
                tn = Get_tn((*h));
                if(read_g->seq[tn].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
                }
                if(read_g->seq[tn].del == 1) continue;
                if(r_flag[tn])
                {
                    C_bases_primary += Get_qe((*h)) - Get_qs((*h));
                }
                else
                {
                    C_bases_alter += Get_qe((*h)) - Get_qs((*h));
                }
            }

            C_bases = C_bases_primary + C_bases_alter;
            if(C_bases_alter < C_bases * ALTER_COV_THRES) continue;

            C_bases = C_bases/R_bases;
            if(C_bases >= asm_opt.recover_atg_cov_min && C_bases <= asm_opt.recover_atg_cov_max)
            {
                available_reads++;
            }
        }


        if(available_reads < (u->n * 0.8) || available_reads == 0)
        {
            u_flag[v] = 0;
        }
        else
        {
            u_flag[v] = 1;
        }

        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            r_flag[rId] = 0;
        }
    }

    free(r_flag);
    return u_flag;
}


void adjust_utg_by_trio(ma_ug_t **ug, asg_t* read_g, uint8_t flag, float drop_rate,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp,
kvec_asg_arc_t_warp* new_rtg_edges)
{
    asg_t* nsg = (*ug)->g;
    uint32_t v, n_vtx = nsg->n_seq;    
    /**
    kvec_t_u32_warp new_rtg_nodes;
    kv_init(new_rtg_nodes.a);
    **/
    purge_dups(*ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, new_rtg_edges, 
        asm_opt.purge_simi_rate, asm_opt.purge_overlap_len, max_hang, min_ovlp, bubble_dist,
        drop_ratio, 1, 1);
    if(asm_opt.recover_atg_cov_min == -1024)
    {
        asm_opt.recover_atg_cov_max = asm_opt.hom_global_coverage/HOM_PEAK_RATE;
        asm_opt.recover_atg_cov_min = asm_opt.recover_atg_cov_max * 0.8;
        ///asm_opt.recover_atg_cov_max = asm_opt.recover_atg_cov_max * 1.2;
        asm_opt.recover_atg_cov_max = INT32_MAX;
    }
    if(asm_opt.recover_atg_cov_max != INT32_MAX)
    {
        fprintf(stderr, "[M::%s] primary contig coverage range: [%d, %d]\n", 
        __func__, asm_opt.recover_atg_cov_min, asm_opt.recover_atg_cov_max);
    }
    else
    {
        fprintf(stderr, "[M::%s] primary contig coverage range: [%d, infinity]\n", 
        __func__, asm_opt.recover_atg_cov_min);
    }
    
    

    ///primary_flag = get_utg_attributes(*ug, read_g, coverage_cut, sources, ruIndex);
    update_unitig_graph((*ug), read_g, coverage_cut, sources, reverse_sources, ruIndex, 0, 
    DOUBLE_CHECK_THRES, flag, drop_rate);

    adjust_utg_advance(read_g, (*ug), reverse_sources, ruIndex);

    nsg = (*ug)->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen((*ug)->u, v) = (*ug)->u.a[v].n;
    }
    clean_trio_untig_graph(*ug, read_g, coverage_cut, sources, reverse_sources, bubble_dist, 
    tipsLen, tip_drop_ratio, stops_threshold, ruIndex, NULL, NULL, 0, 0, 0, 
    chimeric_rate, 0, 0, drop_ratio, flag, drop_rate);

    ///if(flag == MOTHER) fprintf(stderr, "(o.1) c: %u, del: %u, n: %u\n", (*ug)->g->seq[28141].c, (*ug)->g->seq[28141].del, (*ug)->u.a[28141].n);
    
    ///delete_useless_nodes(ug);
    delete_useless_trio_nodes(ug, read_g, coverage_cut, sources, ruIndex);

    update_hap_label(*ug, read_g);

    update_unitig_graph((*ug), read_g, coverage_cut, sources, reverse_sources, ruIndex, 0, 
    DOUBLE_CHECK_THRES, flag, drop_rate);

    renew_utg(ug, read_g, new_rtg_edges);

    if (!(asm_opt.flag & HA_F_BAN_POST_JOIN))
    {
        rescue_missing_overlaps_aggressive(*ug, read_g, sources, coverage_cut, ruIndex, max_hang,
        min_ovlp, 0, 0, 1, NULL);

        renew_utg(ug, read_g, new_rtg_edges);

        rescue_contained_reads_aggressive(*ug, read_g, sources, coverage_cut, ruIndex, max_hang, 
        min_ovlp, 0, 10, 0, 1, NULL, NULL);

        renew_utg(ug, read_g, new_rtg_edges);
    }
    
    update_unitig_graph((*ug), read_g, coverage_cut, sources, reverse_sources, ruIndex, 1, 
    FINAL_DOUBLE_CHECK_THRES, flag, drop_rate);

    update_hap_label(NULL, read_g);

    renew_utg(ug, read_g, new_rtg_edges);

    ///delete_useless_nodes(ug);
    delete_useless_trio_nodes(ug, read_g, coverage_cut, sources, ruIndex);

    if(asm_opt.purge_level_trio == 1)
    {
		purge_dups(*ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, new_rtg_edges, 
        asm_opt.purge_simi_rate, asm_opt.purge_overlap_len, max_hang, min_ovlp, bubble_dist,
        drop_ratio, 1, 0);
		///delete_useless_nodes(ug);
        delete_useless_trio_nodes(ug, read_g, coverage_cut, sources, ruIndex);
	}


    set_drop_trio_flag(*ug);
    
    /**
    kv_destroy(new_rtg_nodes.a);
    **/
}


int debug_untig_length(ma_ug_t *g, uint32_t tipsLen, const char* name)
{
	uint32_t i;

	for (i = 0; i < g->u.n; ++i) {
		ma_utg_t *u = &g->u.a[i];
        if(u->m == 0) continue;
        if(u->n <= tipsLen)
        {
            fprintf(stderr, "i: %u, u->n: %u, tipsLen: %u, %s\n", i, u->n, tipsLen, name);
        }
	}

	return 0;
}



void output_trio_unitig_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
uint8_t flag, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long bubble_dist, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, 
float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp)
{
    char* gfa_name = (char*)malloc(strlen(output_file_name)+100);
    sprintf(gfa_name, "%s.%s.p_ctg.gfa", output_file_name, (flag==FATHER?"hap1":"hap2"));
    fprintf(stderr, "Writing %s to disk... \n", gfa_name);
    FILE* output_file = fopen(gfa_name, "w");

    ma_ug_t *ug = NULL;
    ug = ma_ug_gen(sg);

    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ///print_untig_by_read(ug, "m64011_190830_220126/117834372/ccs", 865264, sources, reverse_sources, "beg");
    adjust_utg_by_trio(&ug, sg, flag, TRIO_THRES, sources, reverse_sources, coverage_cut, 
    bubble_dist, tipsLen, tip_drop_ratio, stops_threshold, ruIndex, chimeric_rate, drop_ratio, 
    max_hang, min_ovlp, &new_rtg_edges);    
    ///debug_utg_graph(ug, sg, 0, 0);
    ///debug_untig_length(ug, tipsLen, gfa_name);
    ///print_untig_by_read(ug, "m64011_190901_095311/125831121/ccs", 2310925, "end");
    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, (flag==FATHER?"h1tg":"h2tg"), output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.%s.p_ctg.noseq.gfa", output_file_name, (flag==FATHER?"hap1":"hap2"));
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, (flag==FATHER?"h1tg":"h2tg"), output_file);
    fclose(output_file);
    if(asm_opt.bed_inconsist_rate != 0)
    {
        sprintf(gfa_name, "%s.%s.p_ctg.lowQ.bed", output_file_name, (flag==FATHER?"hap1":"hap2"));
        output_file = fopen(gfa_name, "w");
        ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
        max_hang, min_ovlp, asm_opt.bed_inconsist_rate, (flag==FATHER?"h1tg":"h2tg"), output_file);
        fclose(output_file);
    }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}


void output_read_graph(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, long long n_read)
{
    fprintf(stderr, "Writing read GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
    sprintf(gfa_name, "%s.read.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_sg_print(sg, &R_INF, coverage_cut, output_file);
    free(gfa_name);
    fclose(output_file);
}


void read_ma(ma_hit_t* x, FILE* fp)
{
    int f_flag;
    f_flag = fread(&(x->qns), sizeof(x->qns), 1, fp);
    f_flag += fread(&(x->qe), sizeof(x->qe), 1, fp);
    f_flag += fread(&(x->tn), sizeof(x->tn), 1, fp);
    f_flag += fread(&(x->ts), sizeof(x->ts), 1, fp);
    f_flag += fread(&(x->te), sizeof(x->te), 1, fp);
    f_flag += fread(&(x->el), sizeof(x->el), 1, fp);
    f_flag += fread(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t;
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->ml = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->rev = t;
    
    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->bl = t;

    f_flag += fread(&(t), sizeof(t), 1, fp);
    x->del = t;
}

int load_ma_hit_ts(ma_hit_t_alloc** x, char* read_file_name)
{
    double startTime = Get_T();
    fprintf(stderr, "Loading ma_hit_ts from disk (%s)... \n", read_file_name);
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        free(index_name);
        return 0;
    }


    long long n_read;
    long long i, k;
    int f_flag;
    f_flag = fread(&n_read, sizeof(n_read), 1, fp);
    (*x) = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*n_read);


    for (i = 0; i < n_read; i++)
    {        
        f_flag += fread(&((*x)[i].is_fully_corrected), sizeof((*x)[i].is_fully_corrected), 1, fp);
        f_flag += fread(&((*x)[i].is_abnormal), sizeof((*x)[i].is_abnormal), 1, fp);
        f_flag += fread(&((*x)[i].length), sizeof((*x)[i].length), 1, fp);
        (*x)[i].size = (*x)[i].length;

        (*x)[i].buffer = NULL;
        if((*x)[i].length == 0) continue;
        
        (*x)[i].buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*(*x)[i].length);

        for (k = 0; k < (*x)[i].length; k++)
        {
            read_ma(&((*x)[i].buffer[k]), fp);
        }  
    }

    free(index_name);    
    fclose(fp);
    fprintf(stderr, "ma_hit_ts has been read (time: %0.2f s).\n", Get_T()-startTime);

    return 1;
}


void write_ma(ma_hit_t* x, FILE* fp)
{
    fwrite(&(x->qns), sizeof(x->qns), 1, fp);
    fwrite(&(x->qe), sizeof(x->qe), 1, fp);
    fwrite(&(x->tn), sizeof(x->tn), 1, fp);
    fwrite(&(x->ts), sizeof(x->ts), 1, fp);
    fwrite(&(x->te), sizeof(x->te), 1, fp);
    fwrite(&(x->el), sizeof(x->el), 1, fp);
    fwrite(&(x->no_l_indel), sizeof(x->no_l_indel), 1, fp);

    uint32_t t = x->ml;
    fwrite(&(t), sizeof(t), 1, fp);
    t = x->rev;
    fwrite(&(t), sizeof(t), 1, fp);

    t = x->bl;
    fwrite(&(t), sizeof(t), 1, fp);
    t =x->del;
    fwrite(&(t), sizeof(t), 1, fp);
}


void write_ma_hit_ts(ma_hit_t_alloc* x, long long n_read, char* read_file_name)
{
    fprintf(stderr, "Writing ma_hit_ts to disk... \n");
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    long long i, k;
    fwrite(&n_read, sizeof(n_read), 1, fp);


    for (i = 0; i < n_read; i++)
    {
        fwrite(&(x[i].is_fully_corrected), sizeof(x[i].is_fully_corrected), 1, fp);
        fwrite(&(x[i].is_abnormal), sizeof(x[i].is_abnormal), 1, fp);
        fwrite(&(x[i].length), sizeof(x[i].length), 1, fp);
        for (k = 0; k < x[i].length; k++)
        {
            write_ma(x[i].buffer + k, fp);
        }  
    }


    free(index_name);
    fflush(fp);    
    fclose(fp);
    fprintf(stderr, "ma_hit_ts has been written.\n");
}

void write_all_data_to_disk(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, All_reads *RNF, char* output_file_name)
{   
	char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
	sprintf(gfa_name, "%s.ec", output_file_name);
	write_All_reads(RNF, gfa_name);

	sprintf(gfa_name, "%s.ovlp.source", output_file_name);
	write_ma_hit_ts(sources, RNF->total_reads, gfa_name);

	sprintf(gfa_name, "%s.ovlp.reverse", output_file_name);
	write_ma_hit_ts(reverse_sources, RNF->total_reads, gfa_name);

	free(gfa_name);
    fprintf(stderr, "bin files have been written.\n");
}

int load_all_data_from_disk(ma_hit_t_alloc **sources, ma_hit_t_alloc **reverse_sources, char* output_file_name)
{
    double startTime = Get_T();
	char* gfa_name = (char*)malloc(strlen(output_file_name)+25);
	sprintf(gfa_name, "%s.ec", output_file_name);
	if (!load_All_reads(&R_INF, gfa_name)) {
		free(gfa_name);
		return 0;
	}
	sprintf(gfa_name, "%s.ovlp.source", output_file_name);
	if (!load_ma_hit_ts(sources, gfa_name)) {
		free(gfa_name);
		return 0;
	}
	sprintf(gfa_name, "%s.ovlp.reverse", output_file_name);
	if (!load_ma_hit_ts(reverse_sources, gfa_name)) {
		free(gfa_name);
		return 0;
	}

    // /////////// hamt containtment bug resolve //////////////
	// sprintf(gfa_name, "%s.ec.was_symm.bin", output_file_name);
    // FILE *fp = fopen(gfa_name, "r");
    // if (!fp){
    //     fprintf(stderr, "symm file missing\n");
    //     free(gfa_name);
    //     return 0;
    // }else{
    //     // paf
    //     for (int i=0; i<R_INF.total_reads; i++){
    //         fread(&R_INF.paf[i].length, sizeof(uint32_t), 1, fp);
    //     }
    //     for (int i=0; i<R_INF.total_reads; i++){
    //         R_INF.paf[i].was_symm = (uint8_t*)calloc(R_INF.paf[i].length, 1);
    //         fread(R_INF.paf[i].was_symm, sizeof(uint8_t), R_INF.paf[i].length, fp);
    //     }
    //     // reverse paf
    //     for (int i=0; i<R_INF.total_reads; i++){
    //         fread(&R_INF.reverse_paf[i].length, sizeof(uint32_t), 1, fp);
    //     }
    //     for (int i=0; i<R_INF.total_reads; i++){
    //         R_INF.reverse_paf[i].was_symm = (uint8_t*)calloc(R_INF.reverse_paf[i].length, 1);
    //         fread(R_INF.reverse_paf[i].was_symm, sizeof(uint8_t), R_INF.reverse_paf[i].length, fp);
    //     }
    //     fclose(fp);
    // }
    

	free(gfa_name);
    fprintf(stderr, "[M::%s] loaded; used %0.2f s\n", __func__, Get_T()-startTime);
	return 1;
}

// count the number of outgoing arcs, excluding reduced arcs
static inline int count_out(const asg_t *g, uint32_t v)
{
	uint32_t i, n, nv = asg_arc_n(g, v);
	const asg_arc_t *av = asg_arc_a(g, v);
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}




// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asg_bub_backtrack_primary(asg_t *g, uint32_t v0, buf_t *b)
{
	uint32_t i, v, qn, tn;
    ///b->S.a[0] is the sink of this bubble
    uint32_t tmp_c = g->seq[b->S.a[0]>>1].c;

	///assert(b->S.n == 1);
	///first remove all nodes in this bubble
	for (i = 0; i < b->b.n; ++i)
    {
        g->seq[b->b.a[i]>>1].c = ALTER_LABLE;
    }

    ///v is the sink of this bubble
	v = b->S.a[0];
	///recover node
	do {
		uint32_t u = b->a[v].p; // u->v
        /****************************may have hap bugs********************************/
		////g->seq[v>>1].c = PRIMARY_LABLE;
        g->seq[v>>1].c = HAP_LABLE;
        /****************************may have hap bugs********************************/
		v = u;
	} while (v != v0);
    ///especially for unitig graph, don't label beg and sink node of a bubble as HAP_LABLE
    ///since in unitig graph, a node may consist of a lot of reads
    g->seq[b->S.a[0]>>1].c = tmp_c;

    ///remove all edges (self/reverse for each edge) in this bubble
	for (i = 0; i < b->e.n; ++i) {
		asg_arc_t *a = &g->arc[b->e.a[i]];
        qn = a->ul>>33;
        tn = a->v>>1;
        if(g->seq[qn].c == ALTER_LABLE && g->seq[tn].c == ALTER_LABLE)
        {
            continue;
        }
		///remove this edge self
		a->del = 1;
		///remove the reverse direction
		asg_arc_del(g, a->v^1, a->ul>>32^1, 1);
	}

    ///v is the sink of this bubble
	v = b->S.a[0];
	///recover node
	do {
		uint32_t u = b->a[v].p; // u->v
		g->seq[v>>1].del = 0;
		asg_arc_del(g, u, v, 0);
		asg_arc_del(g, v^1, u^1, 0);
		v = u;
	} while (v != v0);
}





// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
uint64_t asg_bub_pop1_primary_trio(asg_t *g, ma_ug_t *utg, uint32_t v0, int max_dist, buf_t *b, 
uint32_t positive_flag, uint32_t negative_flag, uint32_t is_pop)
{   
	uint32_t i, n_pending = 0, is_first = 1, cur_m, cur_c, cur_np, cur_nc, to_replace, n_tips, tip_end;
	uint64_t n_pop = 0;
    long long cur_weight = -1, max_weight = -1;
	///if this node has been deleted
	if (g->seq[v0>>1].del || g->seq[v0>>1].c == ALTER_LABLE) return 0; // already deleted
	///if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
    if(get_real_length(g, v0, NULL)<2) return 0;
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = b->a[v0].m = b->a[v0].nc = b->a[v0].np = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);
    n_tips = 0;
    tip_end = (uint32_t)-1;
    uint32_t non_positive_flag = (uint32_t)-1;
    if(positive_flag == FATHER) non_positive_flag = MOTHER;
    if(positive_flag == MOTHER) non_positive_flag = FATHER;

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c, m = b->a[v].m, nc = b->a[v].nc, np = b->a[v].np;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
		///assert(nv > 0);
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
            
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto pop_reset;
            /****************************may have bugs********************************/
            ///important when poping at long untig graph
            if(is_first) l = 0;
            /****************************may have bugs********************************/

			///if this edge has been deleted
			if (av[i].del) continue;

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
                /****************************may have bugs********************************/
                cur_c = cur_m = cur_np = 0; cur_nc = 1;
                if(utg)
                {
                    cur_c = get_num_trio_flag(utg, w>>1, positive_flag);
                    cur_m = get_num_trio_flag(utg, w>>1, negative_flag);
                    cur_np = 0; 
                    if(non_positive_flag != (uint32_t)-1) 
                    {
                        cur_np = get_num_trio_flag(utg, w>>1, non_positive_flag);
                    }
                    cur_nc = utg->u.a[(w>>1)].n;
                }
                

                t->c = c + cur_c;
                t->m = m + cur_m;
                t->nc = nc + cur_nc;
                t->np = np + cur_np;
                /****************************may have bugs********************************/
				///incoming edges of w
				///t->r = count_out(g, w^1);
                t->r = get_real_length(g, w^1, NULL);
				++n_pending;
			} else { // visited before
                /****************************may have bugs********************************/
                cur_c = cur_m = cur_np = 0; cur_nc = 1;
                if(utg)
                {
                    cur_c = get_num_trio_flag(utg, w>>1, positive_flag);
                    cur_m = get_num_trio_flag(utg, w>>1, negative_flag);
                    cur_np = 0; 
                    if(non_positive_flag != (uint32_t)-1) 
                    {
                        cur_np = get_num_trio_flag(utg, w>>1, non_positive_flag);
                    }
                    cur_nc = utg->u.a[(w>>1)].n;
                }
                ///BUG: select the path with less negative_flag, less non_positive_flag, more positive_flag, more distance
                ///FIXED: select the path with less (negative_flag+non_positive_flag), more positive_flag, more distance
                to_replace = 0;

                /****************************may have bugs********************************/
                cur_weight = (long long)(c + cur_c) - ((long long)(m + cur_m) + (long long)(np + cur_np));
                max_weight = (long long)t->c - ((long long)t->m + (long long)t->np);
                if(cur_weight > max_weight)
                {
                    to_replace = 1;
                }
                else if(cur_weight == max_weight)
                {
                    if(nc + cur_nc > t->nc)
                    {
                        to_replace = 1;
                    }
                    else if(nc + cur_nc == t->nc)
                    {
                        if(d + l > t->d)
                        {
                            to_replace = 1;
                        }
                    }
                }
                /****************************may have bugs********************************/

                /**
                if(((m + cur_m) + (np + cur_np)) < (t->m + t->np))
                {
                    to_replace = 1;
                }
                else if(((m + cur_m) + (np + cur_np)) == (t->m + t->np))
                {
                    if(c + cur_c > t->c)
                    {
                        to_replace = 1;
                    }
                    else if(c + cur_c == t->c)
                    {
                        if(nc + cur_nc > t->nc)
                        {
                            to_replace = 1;
                        }
                        else if(nc + cur_nc == t->nc)
                        {
                            if(d + l > t->d)
                            {
                                to_replace = 1;
                            }
                        }
                        
                    }
                }
                **/
                

                if(to_replace)
                {
                    t->p = v;
                    t->m = m + cur_m;
                    t->c = c + cur_c;
                    t->nc = nc + cur_nc;
                    t->np = np + cur_np;
                }
				///c is the weight (is very likely the number of node in this edge) of the parent node
				///select the longest edge (longest meams most reads/longest edge)
                // if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				// if (c + 1 > t->c) t->c = c + 1;
                /****************************may have bugs********************************/
				///update len(v0->w)
                ///node: t->d is not the length from this node's parent
                ///it is the shortest edge
				if (d + l < t->d) t->d = d + l; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
                uint32_t x = get_real_length(g, w, NULL);
                /****************************may have bugs for bubble********************************/
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                /****************************may have bugs for bubble********************************/
				--n_pending;
			}
		}
        is_first = 0;
        //if found a tip
        /****************************may have bugs for bubble********************************/
        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }
        /****************************may have bugs for bubble********************************/
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);


	if(is_pop) asg_bub_backtrack_primary(g, v0, b);

    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = t->m = t->nc = t->np = 0;
	}
	return n_pop;
}


// pop bubbles
int asg_pop_bubble_primary_trio(ma_ug_t *ug, int max_dist, uint32_t positive_flag, uint32_t negative_flag)
{
    asg_t *g = ug->g;
	uint32_t v, n_vtx = g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	//traverse all node with two directions 
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///some node could be deleted
		if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1_primary_trio(ug->g, ug, v, max_dist, &b, positive_flag, negative_flag, 1);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(g);

	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}



// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
uint64_t asg_bub_pop1_primary_trio_debug(ma_ug_t *ug, uint32_t v0, int max_dist, buf_t *b, 
uint32_t positive_flag, uint32_t negative_flag)
{   
    asg_t *g = ug->g;
	uint32_t i, n_pending = 0, is_first = 1, cur_m, cur_c, to_replace, n_tips, tip_end;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del || g->seq[v0>>1].c == ALTER_LABLE) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = b->a[v0].m = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);
    n_tips = 0;
    tip_end = (uint32_t)-1;

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c, m = b->a[v].m;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
		///assert(nv > 0);
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
            
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto pop_reset;
            /****************************may have bugs********************************/
            ///important when poping at long untig graph
            if(is_first) l = 0;
            /****************************may have bugs********************************/

			///if this edge has been deleted
			if (av[i].del) continue;

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l;
                /****************************may have bugs********************************/
                cur_c = get_num_trio_flag(ug, w>>1, positive_flag);
                cur_m = get_num_trio_flag(ug, w>>1, negative_flag);
                t->c = c + cur_c;
                t->m = m + cur_m;
                /****************************may have bugs********************************/
				///incoming edges of w
				t->r = count_out(g, w^1);
				++n_pending;
			} else { // visited before
                /****************************may have bugs********************************/
                cur_c = get_num_trio_flag(ug, w>>1, positive_flag);
                cur_m = get_num_trio_flag(ug, w>>1, negative_flag);
                ///select the way with less negative_flag, more positive_flag, more distance
                to_replace = 0;
                if(m + cur_m < t->m) to_replace = 1;
                if(to_replace == 0 && m + cur_m == t->m && c + cur_c > t->c) to_replace = 1;
                if(to_replace == 0 && m + cur_m == t->m && c + cur_c == t->c && d + l > t->d) to_replace = 1;
                if(to_replace)
                {
                    t->p = v;
                    t->m = m + cur_m;
                    t->c = c + cur_c;
                }
				///c is the weight (is very likely the number of node in this edge) of the parent node
				///select the longest edge (longest meams most reads/longest edge)
                // if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				// if (c + 1 > t->c) t->c = c + 1;
                /****************************may have bugs********************************/
				///update len(v0->w)
                ///node: t->d is not the length from this node's parent
                ///it is the shortest edge
				if (d + l < t->d) t->d = d + l; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
                /****************************may have bugs for bubble********************************/
                /**
                if (x) kv_push(uint32_t, b->S, w);
				///else kv_push(uint32_t, b->T, w); // a tip
                else goto pop_reset;
                **/
                if(x > 0)
                {
                    kv_push(uint32_t, b->S, w);
                }
                else
                {
                    ///at most one tip
                    if(n_tips != 0) goto pop_reset;
                    n_tips++;
                    tip_end = w;
                }
                /****************************may have bugs for bubble********************************/
				--n_pending;
			}
		}
        is_first = 0;
        //if found a tip
        /****************************may have bugs for bubble********************************/
        if(n_tips == 1)
        {
            if(tip_end != (uint32_t)-1 && n_pending == 0 && b->S.n == 0)
            {
                kv_push(uint32_t, b->S, tip_end);
                break;
            }
            else
            {
                goto pop_reset;
            }
        }
        /****************************may have bugs for bubble********************************/
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);

    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = t->m = 0;
	}
	return n_pop;
}



// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t asg_bub_pop1_primary(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0, is_first = 1;
	uint64_t n_pop = 0;
	///if this node has been deleted
	if (g->seq[v0>>1].del || g->seq[v0>>1].c == ALTER_LABLE) return 0; // already deleted
	///asg_arc_n(n0)
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	///S saves nodes with all incoming edges visited
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	///for each node, b->a saves all related information
	b->a[v0].c = b->a[v0].d = 0;
	///b->S is the nodes with all incoming edges visited
	kv_push(uint32_t, b->S, v0);

	do {
		///v is a node that all incoming edges have been visited
		///d is the distance from v0 to v
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///why we have this assert?
		///assert(nv > 0);
		///all out-edges of v
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			/**
			p->ul: |____________31__________|__________1___________|______________32_____________|
								qn            direction of overlap       length of this node (not overlap length)
												(in the view of query)
			p->v : |___________31___________|__________1___________|
								tn             reverse direction of overlap
											(in the view of target)
			p->ol: overlap length
			**/
            
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // v->w with length l
			binfo_t *t = &b->a[w];
			///that means there is a circle, directly terminate the whole bubble poping
			///if (w == v0) goto pop_reset;
            if ((w>>1) == (v0>>1)) goto pop_reset;
            /****************************may have bugs********************************/
            ///important when poping at long untig graph
            if(is_first) l = 0;
            /****************************may have bugs********************************/

			///if this edge has been deleted
			if (av[i].del) continue;

			///push the edge
            ///high 32-bit of g->idx[v] is the start point of v's edges
            //so here is the point of this specfic edge
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			///find a too far path? directly terminate the whole bubble poping
			if (d + l > (uint32_t)max_dist) break; // too far

            ///if this node
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				///t->p is the parent node of 
				///t->s = 1 means w has been visited
				///d is len(v0->v), l is len(v->w), so t->d is len(v0->w)
				t->p = v, t->s = 1, t->d = d + l, t->c = c + 1;
				///incoming edges of w
				t->r = count_out(g, w^1);
				++n_pending;
			} else { // visited before
				///c is the weight (is very likely the number of node in this edge) of the parent node
				///select the longest edge (longest meams most reads/longest edge)
                if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				///update len(v0->w)
                ///node: t->d is not the length from this node's parent
                ///it is the shortest edge
				if (d + l < t->d) t->d = d + l; // update dist
			}
			///assert(t->r > 0);
			//if all incoming edges of w have visited
			//push it to b->S
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				///else kv_push(uint32_t, b->T, w); // a tip
                else goto pop_reset;
				--n_pending;
			}
		}
        is_first = 0;
		///if i < nv, that means (d + l > max_dist)
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);
	asg_bub_backtrack_primary(g, v0, b);
    n_pop = 1;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}


// pop bubbles
int asg_pop_bubble_primary(asg_t *g, int max_dist)
{
	uint32_t v, n_vtx = g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	///set information for each node
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	//traverse all node with two directions 
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		///some node could be deleted
		if (nv < 2 || g->seq[v>>1].del || g->seq[v>>1].c == ALTER_LABLE) continue;
		///some edges could be deleted
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1_primary(g, v, max_dist, &b);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(g);

	if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] popped %lu bubbles\n", __func__, (unsigned long)n_pop);
    }
    return n_pop;
}





int test_triangular_directly(asg_t *g, uint32_t v, 
long long min_edge_length, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    
    uint32_t w;
    int todel;
    long long NodeLen_first[3];
    long long NodeLen_second[3];
    
    asg_arc_t *av = asg_arc_a(g, v);
    if(av[0].v == av[1].v)
    {
        return 0;
    }
    /**********************test first node************************/
    NodeLen_first[0] = NodeLen_first[1] = NodeLen_first[2] = -1;
    if(asg_is_single_edge(g, av[0].v, v>>1) <= 2 && asg_is_single_edge(g, av[1].v, v>>1) <= 2)
    {
        NodeLen_first[asg_is_single_edge(g, av[0].v, v>>1)] = 0;
        NodeLen_first[asg_is_single_edge(g, av[1].v, v>>1)] = 1;
    }
    ///one node has one out-edge, another node has two out-edges
    if(NodeLen_first[1] == -1 || NodeLen_first[2] == -1)
    {
        return 0;
    }
    /**********************test first node************************/
    ///if the potiential edge has already been removed 
    if(av[NodeLen_first[2]].del == 1)
    {
        return 0;
    }

    /**********************test second node************************/
    w = av[NodeLen_first[2]].v^1;
    asg_arc_t *aw = asg_arc_a(g, w);
    uint32_t nw = asg_arc_n(g, w);
    if(nw != 2)
    {
        fprintf(stderr, "[E::%s]error\n", __func__);
    }
    NodeLen_second[0] = NodeLen_second[1] = NodeLen_second[2] = -1;
    if(asg_is_single_edge(g, aw[0].v, w>>1) <= 2 && asg_is_single_edge(g, aw[1].v, w>>1) <= 2)
    {
        NodeLen_second[asg_is_single_edge(g, aw[0].v, w>>1)] = 0;
        NodeLen_second[asg_is_single_edge(g, aw[1].v, w>>1)] = 1;
    }
    ///one node has one out-edge, another node has two out-edges
    if(NodeLen_second[1] == -1 || NodeLen_second[2] == -1)
    {
        return 0;
    }



    todel = 0;
    // if(check_if_diploid(av[0].v, av[1].v, g, reverse_sources, min_edge_length) && 
    // check_if_diploid(aw[0].v, aw[1].v, g, reverse_sources, min_edge_length))
    if(check_if_diploid(av[0].v, av[1].v, g, reverse_sources, min_edge_length, ruIndex) == 1||
    check_if_diploid(aw[0].v, aw[1].v, g, reverse_sources, min_edge_length, ruIndex) == 1)
    {
        todel = 1;
    }
    
    
    
    if(todel)
    {
        ///fprintf(stderr, "v: %u\n", v>>1);
        av[NodeLen_first[2]].del = 1;
        ///remove the reverse direction
        asg_arc_del(g, av[NodeLen_first[2]].v^1, av[NodeLen_first[2]].ul>>32^1, 1);
    }

    return todel;
}



int asg_arc_del_triangular_directly(asg_t *g, long long min_edge_length, 
ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(g, v);
        if (g->seq[v>>1].del)
        {
            continue;
        } 

        if(nv < 2)
        {
            continue;
        }


        n_reduced += test_triangular_directly(g, v, min_edge_length, reverse_sources, ruIndex);
    }

    

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d triangular overlaps, used %.1fs\n", 
                    __func__, n_reduced, Get_T()-startTime);
    return n_reduced;
}




int asg_arc_del_orthology(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;
    uint32_t idx[2];

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);

        // require v to have 2 (still exist) arcs in forward direction
        if (nv < 2 || g->seq[v>>1].del) continue;
        for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
            if (!av[i].del) ++n_arc;
        if (n_arc != 2) continue;

        // collect the 2 arcs
        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                idx[n_arc] = i;
                n_arc++;
            }
        }

        if(check_if_diploid(av[idx[0]].v, av[idx[1]].v, g, reverse_sources, miniedgeLen, ruIndex) == 0)
        {
            float max = av[idx[0]].ol;
            float min = av[idx[1]].ol;

            if(min < drop_ratio * max)
            {
                av[idx[1]].del = 1;
                asg_arc_del(g, av[idx[1]].v^1, av[idx[1]].ul>>32^1, 1);
                n_reduced++;
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    fprintf(stderr, "[M::%s] removed %d different hap overlaps\n", 
    __func__, n_reduced);
    fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);

    return n_reduced;
}




int asg_arc_del_orthology_multiple_way(asg_t *g, ma_hit_t_alloc* reverse_sources, float drop_ratio, 
long long miniedgeLen, R_to_U* ruIndex)
{
    double startTime = Get_T();
    ///the reason is that each read has two direction (query->target, target->query)
    uint32_t v, v_max, v_maxLen, n_vtx = g->n_seq * 2, n_reduced = 0;

    for (v = 0; v < n_vtx; ++v) 
    {
        if (g->seq_vis[v] != 0) continue;
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (nv < 2 || g->seq[v>>1].del) continue;
        n_arc = get_real_length(g, v, NULL);
        if (n_arc < 2) continue;
        v_max = (uint32_t)-1;
        v_maxLen = 0;

        for (i = 0, n_arc = 0; i < nv; i++)
        {
            if (!av[i].del)
            {
                if(v_max == (uint32_t)-1)
                {
                    v_max = av[i].v;
                    v_maxLen = av[i].ol;
                }
                else if(check_if_diploid(v_max, av[i].v, g, reverse_sources, miniedgeLen, ruIndex) == 0)
                {
                    if(av[i].ol < drop_ratio * v_maxLen)
                    {
                        ///fprintf(stderr, "v: %u, v_max: %u, av[%d].v: %u\n", v>>1, v_max>>1, i, av[i].v>>1);
                        
                        av[i].ol = 1;
                        asg_arc_del(g, av[i].v^1, av[i].ul>>32^1, 1);
                        n_reduced++;
                    }
                }
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d different hap overlaps\n", 
        __func__, n_reduced);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    return n_reduced;
}



uint32_t detect_single_path_with_dels_by_length
(asg_t *g, uint32_t begNode, uint32_t* endNode, long long* Len, buf_t* b, long long maxLen)
{
    
    uint32_t v = begNode, w;
    uint32_t kv, kw;
    (*Len) = 0;


    while (1)
    {
        (*Len)++;
        kv = get_real_length(g, v, NULL);
        (*endNode) = v;

        ///if(b) kv_push(uint32_t, b->b, v>>1);
        if(b) kv_push(uint32_t, b->b, v);

        if(kv == 0)
        {
            return END_TIPS;
        }

        if(kv == 2)
        {
            return TWO_OUTPUT;
        }

        if(kv > 2)
        {
            return MUL_OUTPUT;
        }

        if((*Len) > maxLen)
        {
            return LONG_TIPS;
        }

        ///up to here, kv=1
        ///kw must >= 1
        get_real_length(g, v, &w);  // my note: therefore it's safe to use &w (when it's supposed to be an array of at least kv long)
        kw = get_real_length(g, w^1, NULL);
        v = w;
        (*endNode) = v;


        if(kw == 2)
        {
            (*Len)++;
            ///if(b) kv_push(uint32_t, b->b, v>>1);
            if(b) kv_push(uint32_t, b->b, v);
            return TWO_INPUT;
        }

        if(kw > 2)
        {
            (*Len)++;
            ///if(b) kv_push(uint32_t, b->b, v>>1);
            if(b) kv_push(uint32_t, b->b, v);
            return MUL_INPUT;
        }   


        if((v>>1) == (begNode>>1))
        {
            return LOOP;
        } 
        // will stay in this loop & examine more nodes iff the current arc is an 1-to-1 link.
    }

    return LONG_TIPS;   
}



long long asg_arc_del_self_circle_untig(asg_t *g, long long circleLen, int is_drop)
{
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll;
    asg_arc_t *aw;
    uint32_t nw, k;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        if (g->seq[v>>1].del) continue;
        if(is_drop && g->seq[v>>1].c == ALTER_LABLE) continue;

        n_arc = get_real_length(g, v, NULL);
        if (n_arc != 1) continue;

        for (i = 0; i < nv; i++)
        {
            ///actually there is just one un-del edge
            if (!av[i].del)
            {    
                flag = detect_single_path_with_dels_by_length(g, v, &convex, &ll, NULL, circleLen);
                if(ll > circleLen || flag == LONG_TIPS)
                {
                    break;
                }
                if(flag == LOOP)
                {
                    break;
                }
                if(flag != END_TIPS && flag != LONG_TIPS)
                {
                    w = v^1;
                    n_arc = get_real_length(g, w, NULL);
                    if(n_arc == 0)
                    {
                        break;
                    }
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (k = 0; k < nw; k++)
                    {
                        if ((!aw[k].del) && (aw[k].v == (convex^1)))
                        {
                            aw[k].del = 1;
                            asg_arc_del(g, aw[k].v^1, aw[k].ul>>32^1, 1);
                            n_reduced++;
                        }
                    }
                }
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d self-circles\n", 
        __func__, n_reduced);
    }

    return n_reduced;
}


long long get_untig_coverage(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint32_t* b, uint64_t n)
{
    uint64_t k, j;
    uint32_t v;
    ma_hit_t *h;
    long long R_bases = 0, C_bases = 0;
    for (k = 0; k < n; ++k) 
    {
        v = b[k]>>1;
        R_bases += coverage_cut[v].e - coverage_cut[v].s;
        for (j = 0; j < (uint64_t)(sources[v].length); j++)
        {
            h = &(sources[v].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }
    return C_bases/R_bases;
}

/**
void copy_untig(asg_t *g, long long times, uint32_t* b, long long n, C_graph* cg)
{
    if(times <= 1) return;

    times--;
    long long i, j;
    uint64_t tmp;
    for (i = 0; i < times; i++)
    {
        for (j = 0; j < n; j++)
        {
            tmp = 
            kv_push(uint64_t, cg->Node, );
            asg_add_auxiliary_seq_set(g, (b[j]>>1), 0);
        }
    }
}
**/

void init_C_graph(C_graph* g, uint32_t n_seq)
{
    kv_init(g->Nodes);
    kv_init(g->Edges);
    g->pre_n_seq = n_seq;
    g->seqID = n_seq;
}

void destory_C_graph(C_graph* g)
{
    kv_destroy(g->Nodes);
    kv_destroy(g->Edges);
}


long long asg_arc_del_simple_circle_untig(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, asg_t *g, long long circleLen, int is_drop)
{
    // (hamt note) FUNC
    //    treat the following topo for all vertices in g
    //      given v, if v leads to a single path walk that ends at v, cut the arc of v^1->lastSinglePathNode
    //      (v can have more than 2 predecessors)
    // (hamt note) RETURN
    //    total number of arcs cut
    double startTime = Get_T();
    uint32_t v, w, n_vtx = g->n_seq * 2, n_reduced = 0, convex, flag;
    long long ll/**, coverage**/;
    asg_arc_t *aw;
    uint32_t nw, k;
    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    C_graph cg;  // hamt note: this is not used at all
    init_C_graph(&cg, g->n_seq);

    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        ///some node could be deleted
        // (skip if the seq itself deleted / marked)
        if (g->seq[v>>1].del) continue;
        if(is_drop && g->seq[v>>1].c == ALTER_LABLE) continue;
        // (skip if has only 1 or no predecessor)
        if(get_real_length(g, v^1, NULL)<=1) continue;

        n_arc = get_real_length(g, v, NULL);  // number of targets of v
        if (n_arc != 1) continue;

        for (i = 0; i < nv; i++)  // hamt note: this loop is probably just for the look: break is looks nicer than nested if-else
        {
            ///actually there is just one un-del edge
            if (!av[i].del)
            {
                b.b.n = 0;    // reset/init? visited vertices
                flag = detect_single_path_with_dels_by_length(g, v, &convex, &ll, &b, circleLen);  // convex will points to the endNode
                if(ll > circleLen || flag == LONG_TIPS)
                {
                    break;
                }
                if(flag == LOOP)  // the child node of v is in a loop
                {
                    break;
                }
                if(flag != END_TIPS && flag != LONG_TIPS)  // possible tags: TWO/MUL_OUTPUT, TWO/MUL_INPUT
                {
                    // if reach this point, either the endVertex (aka convex) has at least 2 children, or has at least 2 parents
                    if(v == convex && b.b.n > 1)  // single path walk looped back to v, we want to step back 1 vertex
                    {
                        convex = b.b.a[b.b.n - 2];  // -2 because it's index-1, aka n-1-1
                        b.b.n--;
                    }

                    w = v^1;
                    n_arc = get_real_length(g, w, NULL);
                    if(n_arc == 0)
                    {
                        break;
                    }
                    aw = asg_arc_a(g, w);
                    nw = asg_arc_n(g, w);
                    for (k = 0; k < nw; k++)
                    {
                        if ((!aw[k].del) && (aw[k].v == (convex^1)))
                        {
                            // coverage = get_untig_coverage(sources, coverage_cut, b.b.a, b.b.n);
                            // copy_untig(g, (coverage/asm_opt.coverage), b.b.a, b.b.n, &cg);

                            
                            aw[k].del = 1;
                            asg_arc_del(g, aw[k].v^1, aw[k].ul>>32^1, 1);
                            n_reduced++;
                        }
                    }
                }
            }
        }
    }

    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }

    free(b.b.a);
    destory_C_graph(&cg);

    fprintf(stderr, "[M::%s] removed %d self-circles, used %.1fs\n", 
                    __func__, n_reduced, Get_T()-startTime);

    return n_reduced;
}



uint32_t* build_unitig_index(asg_t *sg, ma_ug_t *ug)
{
    if(sg == NULL && ug == NULL) return NULL;
    uint32_t* index = (uint32_t*)calloc(sg->n_seq, sizeof(uint32_t));
    uint32_t i, j;
    for (i = 0; i < ug->u.n; ++i) 
    {
        ma_utg_t *u = &(ug->u.a[i]);
        for (j = 0; j < u->n; j++)
        {
            index[(u->a[j]>>33)] = i;
        }
    }
    return index;
}


int double_check_tangle(uint32_t vBeg, uint32_t vEnd, uint32_t* u_vecs, uint32_t n, asg_t *nsg)
{
    uint32_t v, nv, k, i, j;
    asg_arc_t *av = NULL;
    if(vBeg != (uint32_t)-1)
    {
        v = vBeg;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==(u_vecs[i]>>1)) break;
            }
            ///haven't found
            if(i == n) return 0;
        }
    }

    if(vEnd != (uint32_t)-1)
    {
        v = vEnd^1;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==(u_vecs[i]>>1)) break;
            }
            ///haven't found
            if(i == n) return 0;
        }
    }

    ///scan tangles
    for (j = 0; j < n; j++)
    {
        v = u_vecs[j];
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==(u_vecs[i]>>1)) break;
            }
            if(i == n && av[k].v != (vBeg^1) && av[k].v != vEnd) return 0;
        }



        v = v^1;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==(u_vecs[i]>>1)) break;
            }
            if(i == n && av[k].v != (vBeg^1) && av[k].v != vEnd) return 0;
        }
    }
    
    return 1;
}

int explore_graph(asg_t *nsg, uint32_t vBeg, float single_threshold,
float l_untig_rate_threshold, long long minLongUntig, long long maxShortUntig, uint32_t ignore_d,
buf_t* bb, uint32_t** r, size_t* rm, size_t* rn, uint8_t* visit, ma_ug_t *ug, uint32_t* r_ID)
{
    ///the ID of the end long unitig
    (*r_ID) = (uint32_t)-1;
    uint32_t nv, vBeg_end;
    asg_arc_t *av;

    kvec_t(uint32_t) u_vecs;
    kv_init(u_vecs);
    if(r && rm && rn) kv_reuse(u_vecs, 0, (*rm), (*r));
    
    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    uint32_t vEnd, threshold, num_reads = 0, i, k, v, end = (uint32_t)-1, in = 0, vELen, tmp;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    /**
    bb->b.n = 0;
    threshold = get_long_tip_length(nsg, ut_v, vBeg, &vEnd, bb);
    **/
    bb->b.n = 0;
    if(get_unitig(nsg, ug, vBeg, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb)==LOOP)
    {
        ///the length of LOO is infinite
        kdq_destroy(uint32_t, buf); 
        return 0;
    }
    vBeg_end = vEnd;
    threshold = nodeLen;
    for (i = 0; i < bb->b.n; i++)
    {
        Set_vis(visit, bb->b.a[i], ignore_d);
        Set_vis(visit, bb->b.a[i]^1, ignore_d);
    }    
    v = vBeg;

    kdq_push(uint32_t, buf, v);
    while (kdq_size(buf) != 0)
    {
        in++;
        v = *(kdq_pop(uint32_t, buf));

        ///in == 1 means the start node, it is useless
        if(in != 1)
        {
            ///get current untig length
            /**
            bb->b.n = 0;
            vELen = get_long_tip_length(nsg, ut_v, v, &vEnd, bb);
            **/
            bb->b.n = 0;
            if(get_unitig(nsg, ug, v, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb)==LOOP)
            {
                kdq_destroy(uint32_t, buf);
                return 0;
            }
            vELen = nodeLen;

            ///first long unitig except the start node
            if(if_long_tip_length(nsg, ug, v, &vELen, 
                    minLongUntig, maxShortUntig, l_untig_rate_threshold, threshold)==1)
            {
                if(end == (uint32_t)-1)
                {
                    end = v;
                    continue;
                }
                else
                {
                    in = (uint32_t)-1;
                    break;
                }
            }

            if(vELen > (threshold * single_threshold))
            {
                in = (uint32_t)-1;
                break;
            }

            num_reads = num_reads + vELen;
            if(num_reads > threshold)
            {
                in = (uint32_t)-1;
                break;
            }

            if(r && rm && rn)
            {
                for (i = 0; i < bb->b.n; i++)
                {
                    kv_push(uint32_t, u_vecs, bb->b.a[i]);
                }   
            } 
        }

        tmp = v^1;
        v = vEnd;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;

            if(av[k].v == vBeg)
            {
                in = (uint32_t)-1;
                goto termi;
            }
            
            if(Get_vis(visit,av[k].v,ignore_d)==0)
            {
                kdq_push(uint32_t, buf, av[k].v);

                /**
                bb->b.n = 0;
                get_long_tip_length(nsg, ut_v, av[k].v, &vEnd, bb);
                **/
                bb->b.n = 0;
                get_unitig(nsg, ug, av[k].v, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb);

                for (i = 0; i < bb->b.n; i++)
                {
                    Set_vis(visit, bb->b.a[i], ignore_d);
                }
            }
        }

        ///for start node, we just need one direction
        if(in != 1 && ignore_d)
        {
            v = tmp;
            nv = asg_arc_n(nsg, v);
            av = asg_arc_a(nsg, v);
            for(k = 0; k < nv; k++)
            {
                if(av[k].del) continue;

                if(av[k].v == vBeg)
                {
                    in = (uint32_t)-1;
                    goto termi;
                }

                if(Get_vis(visit,av[k].v,ignore_d)==0)
                {
                    kdq_push(uint32_t, buf, av[k].v);

                    /**
                    bb->b.n = 0;
                    get_long_tip_length(nsg, ut_v, av[k].v, &vEnd, bb);
                    **/
                    bb->b.n = 0;
                    get_unitig(nsg, ug, av[k].v, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb);


                    for (i = 0; i < bb->b.n; i++)
                    {
                        Set_vis(visit, bb->b.a[i], ignore_d);
                    }
                }
            }
        }
    }

    termi:
    kdq_destroy(uint32_t, buf);
    if(r && rm && rn)
    {
        (*rn) = u_vecs.n;
        (*rm) = u_vecs.m;
        (*r) = u_vecs.a;
    }


    (*r_ID) = end;
    ///in == 1 means the end subgraph is the long untig itself
    ///so here is no tangles
    if(in == (uint32_t)-1 || in == 1)
    {
        return 0;
    }
    else
    {
        if(r && rm && rn) return double_check_tangle(vBeg_end, (*r_ID), u_vecs.a, u_vecs.n, nsg);

        return 1;
    }
}


int explore_graph_back(asg_t *nsg, uint32_t start, uint32_t threshold, float single_threshold,
float l_untig_rate_threshold, long long minLongUntig, long long maxShortUntig, uint32_t ignore_d,
uint32_t** r, size_t* rm, size_t* rn, uint8_t* visit, ma_utg_v* ut_v, uint32_t* r_ID)
{
    (*r_ID) = (uint32_t)-1;
    uint32_t nv;
    asg_arc_t *av;

    kvec_t(uint32_t) u_vecs;
    kv_init(u_vecs);
    if(r && rm && rn) kv_reuse(u_vecs, 0, (*rm), (*r));
    
    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    uint32_t num_reads = 0;
    uint32_t v = start, k;
    uint32_t end = (uint32_t)-1;
    uint32_t in = 0;
    Set_vis(visit, v, ignore_d);
    Set_vis(visit, v^1, ignore_d);

    if(nsg->seq[v>>1].del)
    {
        in = (uint32_t)-1;
        goto termi;
    }

    kdq_push(uint32_t, buf, v);
    while (kdq_size(buf) != 0)
    {
        in++;
        v = *(kdq_pop(uint32_t, buf));

        ///in == 1 means the start node, it is useless
        if(in != 1)
        {
            ///first long unitig except the start node
            if(check_long_tip(*ut_v, v>>1, minLongUntig, 
                maxShortUntig, l_untig_rate_threshold, threshold))
            {
                if(end == (uint32_t)-1)
                {
                    end = v;
                    continue;
                }
                else
                {
                    in = (uint32_t)-1;
                    break;
                }
            }

            if(EvaluateLen(*ut_v, v>>1) > (threshold * single_threshold))
            {
                in = (uint32_t)-1;
                break;
            }

            num_reads = num_reads + EvaluateLen(*ut_v, v>>1);
            if(num_reads > threshold)
            {
                in = (uint32_t)-1;
                break;
            }

            if(r && rm && rn) kv_push(uint32_t, u_vecs, v);
        }

        
        
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;

            if(av[k].v == start)
            {
                in = (uint32_t)-1;
                goto termi;
            }
            
            if(Get_vis(visit,av[k].v,ignore_d)==0)
            {
                Set_vis(visit,av[k].v,ignore_d);
                kdq_push(uint32_t, buf, av[k].v);
            }
        }

        ///for start node, we just need one direction
        if(in != 1 && ignore_d)
        {
            v = v^1;
            nv = asg_arc_n(nsg, v);
            av = asg_arc_a(nsg, v);
            for(k = 0; k < nv; k++)
            {
                if(av[k].del) continue;

                if(av[k].v == start)
                {
                    in = (uint32_t)-1;
                    goto termi;
                }

                if(Get_vis(visit,av[k].v,ignore_d)==0)
                {
                    Set_vis(visit,av[k].v,ignore_d);
                    kdq_push(uint32_t, buf, av[k].v);
                }
            }
        }
    }

    termi:
    kdq_destroy(uint32_t, buf);
    if(r && rm && rn)
    {
        (*rn) = u_vecs.n;
        (*rm) = u_vecs.m;
        (*r) = u_vecs.a;
    }
    

    (*r_ID) = end;

    if(in == (uint32_t)-1)
    {
        return 0;
    }
    else
    {
        return 1;
    }
    
}



void output_tangles(uint32_t startID, uint32_t endId, uint32_t* a, uint32_t n, const char* lable)
{
    kvec_t(uint32_t) u_vecs;
    u_vecs.a = a, u_vecs.n = n;

    fprintf(stderr, "\n%sstartID: %u, dir: %u\n", lable, startID>>1, startID&1);
    fprintf(stderr, "%sendID: %u, dir: %u\n", lable, endId>>1, endId&1);
    uint32_t ijk;
    for (ijk = 0; ijk < u_vecs.n; ijk++)
    {
        fprintf(stderr, "%stangleID: %u, dir: %u\n", lable, u_vecs.a[ijk]>>1, u_vecs.a[ijk]&1);
    }
}


/**
 * @par result Is the address of an allocated struct. It gets a copy. 
 *             Use hamt_get_arc to return pointer to the original asg_arc_t item.
 * @return 1 if found, 0 not found
*/
int get_arc(asg_t *g, uint32_t src, uint32_t dest, asg_arc_t* result)
{
    uint32_t i;
    if(g->seq[src>>1].del) return 0;

    uint32_t nv = asg_arc_n(g, src);
    asg_arc_t *av = asg_arc_a(g, src);
    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if(av[i].v == dest)
        {
            (*result) = av[i];
            return 1;
        } 
    }

    return 0;
}

uint32_t insert_index(asg_t *g, uint32_t v)
{
    uint32_t v_tx = g->n_seq * 2;
    if(v >= v_tx) return (uint32_t)-1;
    if(asg_arc_n(g, v)!=0) return (g->idx[v]>>32) + asg_arc_n(g, v);
    ///now v itself does not have any edge
    while (v < v_tx && asg_arc_n(g, v) == 0){v++;}
    ///means there are no edge at the whole graph
    if(v>=v_tx) return g->n_arc;
    return (g->idx[v]>>32);
}

asg_arc_t* insert_index_p(asg_t *g, long long index, uint32_t m_distance)
{
    ///each edge has two direction
    if (g->n_arc + m_distance > g->m_arc) 
    {
        ///g->m_arc = g->n_arc + (m_distance<<1);
        g->m_arc = (g->n_arc + m_distance)<<1;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
	}
    long long i = g->n_arc; i--;
    for (; i >= index; i--)
    {
        g->arc[i+m_distance] = g->arc[i];
    }

    g->n_arc = g->n_arc + m_distance;

    return &(g->arc[index]);
}

void exchange_arcs(asg_arc_t* x, asg_arc_t* y)
{
    asg_arc_t k;
    k = (*x);(*x) = (*y);(*y) = k;
}

void insert_arc(asg_t *g, long long index, uint32_t m_distance, uint32_t src, uint32_t srcLen, uint32_t dest, 
uint32_t oLen, uint8_t strong, uint8_t el, uint8_t no_l_indel)
{
    asg_arc_t* p;
    uint32_t l;
    long long i;
    p = insert_index_p(g, index, m_distance);
    p->del = !!(0);p->el = el; p->no_l_indel = no_l_indel; p->strong = strong; p->ol = oLen;
    p->v = dest;
    p->ul = src; p->ul = p->ul << 32; l = srcLen - oLen; p->ul = p->ul | l;
    for (i = index - 1; i >= 0 && p->ul <= g->arc[i].ul; i--)
    {
        if((!g->arc[i].del)&&(g->arc[i].v==p->v) && ((g->arc[i].ul>>32)==(p->ul>>32)))
        {
            p->del = !!(1);
            break;
        }
        exchange_arcs(p, &(g->arc[i]));p = &(g->arc[i]);
    }
    

    if(asg_arc_n(g, src) == 0)
    {
        g->idx[src] = index; g->idx[src] = g->idx[src]<<32; g->idx[src] = g->idx[src] | 1;
    }
    else
    {
        g->idx[src] = g->idx[src] + 1;
    }
      
    uint32_t v, v_tx = g->n_seq*2;
    uint64_t add = 1; add = add << 32;
    for (v = src+1; v < v_tx; v++)
    {
        if(asg_arc_n(g, v) == 0) continue;
        g->idx[v] += add;
    }
}

int asg_append_edges_to_srt(asg_t *g, uint32_t src, 
uint32_t srcLen, uint32_t dest, uint32_t oLen,
uint8_t strong, uint8_t el, uint8_t no_l_indel)
{
    uint32_t v_tx = g->n_seq * 2;
    /**
    uint32_t d_i = 0, current_i, next_i, src_n, s_index, e_index;
    while(d_i < v_tx)
    {
        current_i = d_i;
        d_i++;
        src_n = asg_arc_n(g, current_i);
        if(src_n == 0) continue;
        s_index = (g->idx[current_i]>>32);
        e_index = s_index + src_n;

        while (d_i < v_tx && asg_arc_n(g, d_i) == 0){d_i++;}
        if(d_i>=v_tx) break;
        next_i = d_i;


        if(current_i != v_tx)
        {
            if(e_index != (g->idx[next_i]>>32))
            {
                fprintf(stderr, "v_tx: %u, current_i: %u, node: %u, s_index: %u, e_index: %u, g->idx[next_i]>>32: %u\n",
                v_tx, current_i, current_i>>1, s_index, e_index, g->idx[next_i]>>32);
            }
        }
        else 
        {
            if(e_index != g->n_arc)
            {
                fprintf(stderr, "ERROR2\n");
            }
        }
    }

    for (d_i = 0; d_i < v_tx; d_i++)
    {
        uint32_t nv = asg_arc_n(g, d_i), s_i;
        asg_arc_t *av = asg_arc_a(g, d_i);
        asg_arc_t forward, backward;
        for (s_i = 0; s_i < nv; s_i++)
        {
            if(av[s_i].del) continue;
            forward = av[s_i];
            if(get_arc(g, forward.ul>>32, forward.v, &forward) == 0)
            {
                fprintf(stderr, "ERROR1\n");
            }
            if(forward.del != av[s_i].del || forward.el != av[s_i].el || 
               forward.no_l_indel != av[s_i].no_l_indel ||forward.ol != av[s_i].ol || 
               forward.strong != av[s_i].strong || forward.ul != av[s_i].ul || 
               forward.v != av[s_i].v)
            {
                fprintf(stderr, "ERROR2\n");
            }

            if(get_arc(g, forward.v^1, (forward.ul>>32)^1, &backward) == 0)
            {
                fprintf(stderr, "ERROR3\n");
            }

            if(forward.ol != backward.ol)
            {
                fprintf(stderr, "forward.ol: %u, backward.ol: %u\n", 
                forward.ol, backward.ol);
            }

        }
    }
    **/


    if (src >= v_tx || dest >= v_tx)
    {
        return 0;
    }

    ///we need to link src--->dest and (dest^1)----->(src^1)
    uint32_t src_i = insert_index(g, src);
    insert_arc(g, src_i, 1, src, srcLen, dest, oLen, strong, el, no_l_indel);


    /**
    if (src >= v_tx || (src^1) >= v_tx || dest >= v_tx || (dest^1) >= v_tx)
    {
        return 0;
    }

    ///we need to link src--->dest and (dest^1)----->(src^1)
    uint32_t src_i = insert_index(g, src);
    insert_arc(g, src_i, 1, src, srcLen, dest, oLen);
    uint32_t dest_i = insert_index(g, (dest^1));
    insert_arc(g, dest_i, 1, dest^1, destLen, src^1, oLen);
    **/
    
    
    ///we need to link src--->dest and (dest^1)----->(src^1)
    /**
    if(src > (dest^1))
    {
        uint32_t src_i = insert_index(g, src);
        insert_arc(g, src_i, src, srcLen, dest, oLen);
        uint32_t dest_i = insert_index(g, (dest^1));
        insert_arc(g, dest_i, dest^1, destLen, src^1, oLen);
    }
    else
    {
        uint32_t dest_i = insert_index(g, (dest^1));
        insert_arc(g, dest_i, dest^1, destLen, src^1, oLen);
        uint32_t src_i = insert_index(g, src);
        insert_arc(g, src_i, src, srcLen, dest, oLen);
    }
    **/

   return 1;
}

void append_ma_utg_t(ma_utg_t* v_x, ma_utg_t* v_y)
{
    if(v_x->m < (v_x->n + v_y->n))
    {
        v_x->m = v_x->n + v_y->n;
        v_x->a = (uint64_t*)realloc(v_x->a, v_x->m * sizeof(uint64_t));
    }
    memcpy(v_x->a+v_x->n, v_y->a, v_y->n*sizeof(uint64_t));
    v_x->n = v_x->n + v_y->n;
    free(v_y->a);
    v_y->m=v_y->n=0;v_y->a=NULL;
}


#define UNROLL 0
#define CONVEX 1
void merge_nodes(ma_ug_t *ug, uint32_t startID, uint32_t endId, uint32_t* a, uint32_t n, 
uint32_t type)
{
    ma_utg_t *v_x = NULL, *v_y = NULL;
    asg_t* nsg = ug->g;
    kvec_t(uint32_t) u_vecs;
    uint32_t i = 0, maxEvaluateLen = 0, maxBaseLen = 0, v, totalEvaluateLen = 0;
    u_vecs.a = a, u_vecs.n = n;
    if(u_vecs.n > 0) 
    {
        v_x = &(ug->u.a[u_vecs.a[0]>>1]);
        asg_seq_del(nsg, u_vecs.a[0]>>1);
        maxEvaluateLen = EvaluateLen(ug->u, u_vecs.a[0]>>1);
        maxBaseLen = v_x->len;
        totalEvaluateLen += EvaluateLen(ug->u, u_vecs.a[0]>>1);
    }

    
    for (i = 1; i < u_vecs.n; i++)
    {
        v_y = &(ug->u.a[u_vecs.a[i]>>1]);
        if(maxEvaluateLen <= EvaluateLen(ug->u, u_vecs.a[i]>>1))
        {
            maxEvaluateLen = EvaluateLen(ug->u, u_vecs.a[i]>>1);
            maxBaseLen = v_y->len;
        }
        totalEvaluateLen += EvaluateLen(ug->u, u_vecs.a[i]>>1);

        
        append_ma_utg_t(v_x, v_y);
        asg_seq_del(nsg, u_vecs.a[i]>>1);
    }


    if(u_vecs.n > 0)
    {
        i = 0;
        nsg->seq[u_vecs.a[0]>>1].del = !!(0);


        EvaluateLen(ug->u, u_vecs.a[0]>>1) = maxEvaluateLen;
        totalEvaluateLen = totalEvaluateLen / 2;
        if(totalEvaluateLen > EvaluateLen(ug->u, u_vecs.a[0]>>1))
        {
            EvaluateLen(ug->u, u_vecs.a[0]>>1) = totalEvaluateLen;
        }
        IsMerge(ug->u, u_vecs.a[0]>>1)++;


        v_x->len = maxBaseLen;
        /****************************may have bugs********************************/
        nsg->seq[u_vecs.a[0]>>1].len = maxBaseLen;
        /****************************may have bugs********************************/
        v = u_vecs.a[0];

        if(type == UNROLL)
        {
            if(startID != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, startID, ug->u.a[startID>>1].len, v, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, v^1, ug->u.a[v>>1].len, startID^1, 0, 0, 0, 0);
            }
            
            if(endId != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, endId, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, endId^1, ug->u.a[endId>>1].len, v^1, 0, 0, 0 ,0);
            }
        }

        if(type == CONVEX)
        {
            if(startID != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, startID^1, 0, 0, 0, 0);
                asg_append_edges_to_srt(nsg, startID, ug->u.a[startID>>1].len, v^1, 0, 0, 0 ,0);
            }

            if(endId != (uint32_t)-1)
            {
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, endId, 0, 0, 0 ,0);
                asg_append_edges_to_srt(nsg, endId^1, ug->u.a[endId>>1].len, v^1, 0, 0, 0 ,0);
            }
        }
    }
}

///return 1 is what we want
///as for return value: 0: do nothing, 1: unroll, 2: convex
#define CONVEX_M 1
#define UNROLL_M 2
#define UNROLL_E 3
inline uint32_t walk_through(asg_t *read_g, ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, buf_t* b_0, buf_t* b_1,
kvec_t_u32_warp* u_vecs, uint8_t* visit, uint32_t v, uint32_t* r_beg, uint32_t* r_end, 
uint32_t* r_next_uID, R_to_U* ruIndex)
{
    (*r_beg) = (*r_end) = (uint32_t)-1;
    asg_t* nsg = ug->g;
    uint32_t i, beg, end, primaryLen, returnFlag;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;

    /*****************************simple checking**********************************/
    beg = v;
    if(nsg->seq[beg>>1].del || asg_arc_n(nsg, beg) <= 0 || get_real_length(nsg, beg, NULL)<=0)
    {
        return 0;
    }
    //  hamt note:let beg not be in the middle of a single path,
    //            or the beginning of that single path has multiple targets
    if(get_real_length(nsg, beg^1, NULL) == 1)
    {
        get_real_length(nsg, beg^1, &end);
        if(get_real_length(nsg, end^1, NULL) == 1)
        {
            return 0;
        }
    }

    ///if the contig here is too small
    /**
    primaryLen = get_long_tip_length(nsg, &(ug->u), beg, &end, NULL);
    if(primaryLen == 0 || primaryLen < minLongUntig || ug->u.a[beg>>1].circ)
    {
        return 0;
    }
    **/
    primaryLen = get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, NULL);
    if(primaryLen == LOOP || nodeLen < minLongUntig)
    {
        return 0;
    } 
    primaryLen = nodeLen;
    
    if(get_real_length(nsg, end, NULL) <= 0)
    {
        return 0;
    }
    /*****************************simple checking**********************************/
    (*r_beg) = beg; (*r_end) = end;
    
    /*****************************adjacent checking**********************************/
    uint32_t nw = asg_arc_n(nsg, end), n_arc = 0, next_uID, next_uID_verify;
    asg_arc_t *aw = asg_arc_a(nsg, end);
    for (i = 0; i < nw; i++)
    {
        if(!aw[i].del)
        {
            n_arc++;
            ////we don't want any long tip here
            if(if_long_tip_length(nsg, ug, aw[i].v, NULL, 
                minLongUntig, maxShortUntig, l_untig_rate, primaryLen)==1)
            {
                n_arc = 0;
                break;
            }
        }
    }

    if(n_arc == 0)
    {
        return 0;
    }
    /*****************************adjacent checking**********************************/

    ///here all out-nodes of v are small untigs
    if(explore_graph(nsg, beg, max_node_threshold, l_untig_rate, 
    minLongUntig, maxShortUntig, 1, b_0, &(u_vecs->a.a), &(u_vecs->a.m), 
    &(u_vecs->a.n), visit, ug, &next_uID) == 1)
    {
        (*r_next_uID) = next_uID;
        if(next_uID != (uint32_t)-1 && u_vecs->a.n == 1 
                        && IsMerge(ug->u, u_vecs->a.a[0]>>1) > 0)
        {
            return 0;
        }
        ///if next_uID == (uint32_t)-1, that means we found an end subgraph
        if(next_uID != (uint32_t)-1)
        {
            ///need to avoid missassembly
            ///merge all tangle as two types: 1. convex; 2, unroll them as a line
            ///should do somthing here, but we just skip it for covenience
            returnFlag = explore_graph(nsg, beg, max_node_threshold, l_untig_rate, 
            minLongUntig, maxShortUntig, 0, b_0, NULL, NULL, NULL, visit, ug, 
            &next_uID_verify);

            if((returnFlag == 1 && next_uID_verify != next_uID) || (returnFlag == 0))
            {
                //output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)("###"));
                returnFlag = 0;
            }
            else
            {
                returnFlag = 1;
            }
            
            // #define UNAVAILABLE (uint32_t)-1
            // #define PLOID 0
            // #define NON_PLOID 1
            /**
            if(returnFlag == 1 && check_if_diploid_untigs(nsg, read_g, beg, next_uID, &(ug->u),
                reverse_sources, minLongUntig-1, 0.3, b_0, b_1, 0, ruIndex) == 1)**/
            if(returnFlag == 1 && check_different_haps(nsg, ug, read_g, beg, next_uID, 
                reverse_sources, b_0, b_1, ruIndex, minLongUntig-1, 1) == PLOID)
            {
                ///output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)("???"));
                returnFlag = 0;
            }

            if(returnFlag == 0)
            {
                merge_nodes(ug, end, next_uID, u_vecs->a.a, u_vecs->a.n, CONVEX);
                return CONVEX_M;
            }
        }

        ///actually it is not possible here
        if(u_vecs->a.n <= 0 || next_uID>>1 == beg>>1)
        {
            return 0;
        }
        ///output_tangles(beg, next_uID, u_vecs->a.a, u_vecs->a.n, (char*)(""));
        merge_nodes(ug, end, next_uID, u_vecs->a.a, u_vecs->a.n, UNROLL);
        if(next_uID != (uint32_t)-1)
        {
            (*r_next_uID) = next_uID;
            return UNROLL_M;
        }
        else ///end tangle
        {
            if(u_vecs->a.n > 0)
            {
                (*r_next_uID) = u_vecs->a.a[0];
            } 
            return UNROLL_E;
        }
    }

    return 0;
}


void adjust_asg_by_ug(ma_ug_t *ug, asg_t *read_g)
{
    uint32_t v, n_vtx, i = 0, qn;
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    ma_utg_t* node = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == ALTER_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                read_g->seq[qn].c = ALTER_LABLE;
            }
        }
    }


    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == ALTER_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                ///read_g->seq[qn].c = ALTER_LABLE;
                asg_seq_drop(read_g, qn);
            }
        }
    }

    asg_cleanup(read_g);
	asg_symm(read_g);
}

void lable_hap_asg_by_ug(ma_ug_t *ug, asg_t *read_g)
{
    uint32_t v, n_vtx, i = 0, qn;
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    ma_utg_t* node = NULL;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].c == HAP_LABLE)
        {
            node = &(ug->u.a[v]);
            for (i = 0; i < node->n; i++)
            {
                qn = node->a[i]>>33;
                read_g->seq[qn].c = HAP_LABLE;
            }
        }
    }

}

///qn is the read Id
///self_offset is the offset of this read in contig
void query_reverse_sources(asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t qn, uint32_t self_offset, kvec_t_u64_warp* u_vecs)
{
    uint32_t i, rId, is_Unitig, cId;
    uint64_t mode;
    ///means the reads coming from different haplotype have already been purged
    if(read_g->seq[qn].c == HAP_LABLE)
    {
        cId = (uint32_t)-1;
        mode = self_offset; mode = mode << 33; mode = mode|cId; mode = mode | (uint64_t)(0x100000000);
        kv_push(uint64_t, u_vecs->a, mode);
        return;
    }


    for (i = 0; i < reverse_sources[qn].length; i++)
    {
        rId = Get_tn(reverse_sources[qn].buffer[i]);
        ///there are three cases: 
        ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
        ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
        ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
        if(read_g->seq[rId].del == 1)
        {
            ///get the id of read that contains it 
            get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
            if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
        }

        ///here rId is the id of the read coming from the different haplotype
        ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
        get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
        if(is_Unitig == 0) continue;
        ///if the read is at alternative contigs, cId might be (uint32_t)-1
        ///if(cId == (uint32_t)-1)
        mode = self_offset; mode = mode << 33; mode = mode|(uint64_t)(cId);
        kv_push(uint64_t, u_vecs->a, mode);
    }
    
}



uint32_t get_rId_from_contig_by_offset(ma_ug_t *ug, rIdContig* array, uint32_t offset)
{
    uint32_t uId, rId;
    ma_utg_t* reads;
    for (;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (;array->readI < reads->n; array->readI++, array->offset++)
        {
            if(array->offset == offset)
            {
                rId = reads->a[array->readI]>>33;
                return rId;
            }
        }
        array->readI = 0;
    }


    array->offset = 0;
    for (array->untigI = 0;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (array->readI = 0;array->readI < reads->n; array->readI++, array->offset++)
        {
            if(array->offset == offset)
            {
                rId = reads->a[array->readI]>>33;
                return rId;
            }
        }
    }

    array->untigI = array->readI = array->offset = 0;
    return (uint32_t)-1;
}

inline uint32_t get_contig_len(ma_ug_t *ug, buf_t* b_0)
{
    uint32_t uId, untigI, Len = 0;
    ma_utg_t* reads;
    
    for (untigI = 0;untigI < b_0->b.n; untigI++)
    {
        uId = b_0->b.a[untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        Len = Len + reads->n;
    }

    return Len;
}

uint32_t get_offset_from_contig_by_rId(ma_ug_t *ug, rIdContig* array, uint32_t query_rId)
{
    uint32_t uId, rId;
    ma_utg_t* reads;
    for (;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (;array->readI < reads->n; array->readI++, array->offset++)
        {
            rId = reads->a[array->readI]>>33;
            if(rId == query_rId) return array->offset;
        }
        array->readI = 0;
    }


    array->offset = 0;
    for (array->untigI = 0;array->untigI < array->b_0->b.n; array->untigI++)
    {
        uId = array->b_0->b.a[array->untigI]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (array->readI = 0;array->readI < reads->n; array->readI++, array->offset++)
        {
            rId = reads->a[array->readI]>>33;
            if(rId == query_rId) return array->offset;
        }
    }

    array->untigI = array->readI = array->offset = 0;
    return (uint32_t)-1;
}

///tn is the cId, tn_off is the order of rId with the same cId
uint32_t get_reverseId(asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, uint32_t qn, uint32_t tn, uint32_t tn_off)
{
    uint32_t i, rId, is_Unitig, cId, cId_off = 0;
    for (i = 0; i < reverse_sources[qn].length; i++)
    {
        rId = Get_tn(reverse_sources[qn].buffer[i]);
        ///there are three cases: 
        ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
        ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
        ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
        if(read_g->seq[rId].del == 1)
        {
            ///get the id of read that contains it 
            get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
            if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
        }

        ///here rId is the id of the read coming from the different haplotype
        ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
        get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
        if(is_Unitig == 0) continue;
        ///if the read is at alternative contigs, cId might be (uint32_t)-1
        ///if(cId == (uint32_t)-1)
        if(cId == tn)
        {
            if(cId_off == tn_off) return rId;
            cId_off++;
        }
    }

    return (uint32_t)-1;
}

uint32_t get_contig_overlap_interval(uint32_t is_reverse, 
uint32_t q_beg, uint32_t q_end, uint32_t qLen, uint32_t t_beg, uint32_t t_end, uint32_t tLen, 
uint32_t* r_q_beg, uint32_t* r_q_end, uint32_t* r_t_beg, uint32_t* r_t_end)
{
    uint32_t k;
    (*r_q_beg) = (*r_q_end) = (*r_t_beg) = (*r_t_end) = (uint32_t)-1;
    if(q_end >= qLen || t_end >= tLen) return 0;
    if(q_beg >= qLen || t_beg >= tLen) return 0;

    if(is_reverse)
    {
        ///k = q_beg; q_beg = q_end; q_end = k;
        ///q_beg = qLen - q_beg - 1; q_end = qLen - q_end - 1; k = q_beg; q_beg = q_end; q_end = k;
        ///k = t_beg; t_beg = t_end; t_end = k;
        t_beg = tLen - t_beg - 1; t_end = tLen - t_end - 1; k = t_beg; t_beg = t_end; t_end = k;
    }




    if(q_beg <= t_beg)
    {
        t_beg = t_beg - q_beg; q_beg = 0;
    }
    else
    {
        q_beg = q_beg - t_beg; t_beg = 0;
    }

    uint32_t q_right_length = qLen - q_end - 1;
    uint32_t t_right_length = tLen - t_end - 1;
    if(q_right_length <= t_right_length)
    {
        q_end = qLen - 1; t_end = t_end + q_right_length;
    }
    else
    {
        t_end = tLen - 1; q_end = q_end + t_right_length;
    }

    if(is_reverse)
    {
        ///q_beg = qLen - q_beg - 1; q_end = qLen - q_end - 1; k = q_beg; q_beg = q_end; q_end = k;
        t_beg = tLen - t_beg - 1; t_end = tLen - t_end - 1; k = t_beg; t_beg = t_end; t_end = k;
    }

    (*r_q_beg) = q_beg; (*r_q_end) = q_end;
    (*r_t_beg) = t_beg; (*r_t_end) = t_end;
    return 1;
}


int get_haplotype_rate(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, buf_t* b_0, uint32_t beg, uint32_t end, uint32_t query_cId, float match_rate)
{
    uint32_t i, j, k, num, match_num, offset, uId, qn, rId, is_Unitig, cId;
    ma_utg_t* reads;

    offset = match_num = num = 0;
    for (i = 0; i < b_0->b.n; i++)
    {
        uId = b_0->b.a[i]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        for (j = 0; j < reads->n; j++, offset++)
        {

            if(offset < beg || offset > end) continue;
            qn = reads->a[j]>>33;
            if(reverse_sources[qn].length > 0) num++;
            for (k = 0; k < reverse_sources[qn].length; k++)
            {
                rId = Get_tn(reverse_sources[qn].buffer[k]);
                ///there are three cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                ///3. read has bee deleted, get_R_to_U() return the id of read that contains it 
                if(read_g->seq[rId].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                    if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
                }

                ///here rId is the id of the read coming from the different haplotype
                ///cId is the id of the corresponding contig (note here is the contig, instead of untig)
                get_R_to_U(ruIndex, rId, &cId, &is_Unitig);
                if(is_Unitig == 0) continue;
                ///if the read is at alternative contigs, cId might be (uint32_t)-1
                
                if(cId == query_cId)
                {
                    match_num++;
                    break;
                } 
                
            }
        }
    } 

    if(match_num >= num*match_rate) return 1;
    return 0; 
}

#define fully_cover_rate 0.7
#define extraord_rate 0.2
///#define hap_seed 20
#define hap_seed 5
#define GetCid(x) ((uint64_t)(0x7fffffff)&(uint64_t)(x))
#define GetOff(x) ((uint64_t)(x)>>33)
#define IfAlter(x) (GetCid((x))==(0x7fffffff))
#define IfColor(x) (((uint64_t)(0x100000000)&(uint64_t)(x))!=0)
#define IfVisit(x) (((uint64_t)(0x80000000)&(uint64_t)(x))!=0)
#define SetVisit(x) ((x) = ((uint64_t)(0x80000000)|(uint64_t)(x)))
inline int get_useful_contig_advance(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, kvec_t_u64_warp* r_vecs, buf_t* b_0, uint64_t* contigBeg, asg_t *bi_g, 
uint32_t currentId, float Hap_rate, uint32_t long_hap_overlap, float long_hap_overlap_rate,
uint32_t is_bi_edge)
{
    uint32_t i = 0, j, k, sLen, is_found, rId, beg, end, is_build;
    uint32_t q_beg, q_end, t_beg, t_end, qLen, tLen;
    uint32_t interval_q_beg, interval_q_end, interval_t_beg, interval_t_end;
    uint64_t anchor_offset, offset;
    uint64_t anchor_cId, cId;

    buf_t b_tid;
    memset(&b_tid, 0, sizeof(buf_t));

    kvec_t(Hap_Align) seed;
    kv_init(seed);

    kvec_t(Hap_Align) alignment;
    kv_init(alignment);

    Hap_Align x;
    rIdContig iterator, iterator_tid;
    iterator.b_0 = b_0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    if(r_vecs->a.n == 0) return -1;
    ///number of reads in this contig
    qLen = get_contig_len(ug, b_0);

    i = 0;
    alignment.n = 0;
    while (i < r_vecs->a.n)
    {
        ///there are three cases:
        /**
         * 1) u_vecs->a.a[i] has haplotype infor at the primary contigs
         * 2) u_vecs->a.a[i] has haplotype infor at the alternative contigs, 
         *    IfAlter(u_vecs->a.a[i])==1
         * 3) u_vecs->a.a[i] itself is labled with HAP_LABLE, 
         *    IfColor(u_vecs->a.a[i])==1
         * 4) if u_vecs->a.a[i] has been labled
         *    IfVisit(u_vecs->a.a[i])==1
         * 4) u_vecs->a.a[i] does not have any haplotype infor. In this case, 
         *    self_offsetLen is not consecutive 
        **/
        if(IfAlter(r_vecs->a.a[i]) || IfColor(r_vecs->a.a[i]) || IfVisit(r_vecs->a.a[i]))
        {
            i++;
            continue;
        }
        
        ///anchor_offset is the self offset, instead of the offset in anchor_cId
        anchor_offset = GetOff(r_vecs->a.a[i]); anchor_cId = GetCid(r_vecs->a.a[i]); 
        sLen = 0; is_found = 0; seed.n = 0;
        for (j = i; j < r_vecs->a.n; j++)
        {
            ///offset is the self offset at query
            ///cId is the target Id (contig)
            offset = GetOff(r_vecs->a.a[j]); 
            cId = GetCid(r_vecs->a.a[j]);

            ///means we go ahead one step at query
            if(offset != anchor_offset)
            {
                ///sLen must be >=1
                if(is_found == 0) 
                {
                    break;
                }
                is_found = 0;
                sLen++;
                anchor_offset = offset;
            } 

            if(IfVisit(r_vecs->a.a[j]) || IfColor(r_vecs->a.a[j])) continue;

            if(cId == anchor_cId)
            {
                x.q_pos = GetOff(r_vecs->a.a[j]);
                x.t_id = GetCid(r_vecs->a.a[j]);
                x.is_color = 1;
                x.t_pos = (uint32_t)-1;
                SetVisit(r_vecs->a.a[j]);
                is_found++;
                kv_push(Hap_Align, seed, x);
            } 
        }

        if(j == r_vecs->a.n && is_found > 0)
        {
            sLen++;
        }


        ///don't need to know how long of this seed
        if(sLen >= hap_seed && seed.n > 0)
        {
            uint32_t inner_off = 0, pre_q_pos = (uint32_t)-1, RrId;
            iterator.offset = iterator.readI = iterator.untigI = 0;

            ///for one vector, all t_id should be same; that is why we recover tid here
            beg = (uint32_t)(contigBeg[seed.a[0].t_id]);
            b_tid.b.n = 0;ug->u.a[beg>>1].circ = 0;
            get_long_tip_length(ug->g, &(ug->u), beg, &end, &b_tid);
            iterator_tid.b_0 = &b_tid;
            iterator_tid.offset = iterator_tid.readI = iterator_tid.untigI = 0;
            ///we have already got qLen at the begining
            tLen = get_contig_len(ug, &b_tid);

            ///just set t_pos
            for (k = 0; k < seed.n; k++)
            {
                rId = get_rId_from_contig_by_offset(ug, &iterator, seed.a[k].q_pos);
                ///actually we shouldn't have this case
                if(rId == (uint32_t)-1) continue;
                if(pre_q_pos != seed.a[k].q_pos)
                {
                    inner_off = 0;
                    pre_q_pos = seed.a[k].q_pos;
                }
                else
                {
                    inner_off++;
                }
                ///rId is the read ID at the
                RrId = get_reverseId(read_g, reverse_sources, ruIndex, rId, seed.a[k].t_id, 
                inner_off);
                ///if(RrId == (uint32_t)-1) fprintf(stderr, "ERROR1\n");
                ///for one vector, all t_id should be same
                seed.a[k].t_pos = get_offset_from_contig_by_rId(ug, &iterator_tid, RrId);
                // if(seed.a[k].t_pos == (uint32_t)-1) fprintf(stderr, "ERROR2\n");
                // uint32_t debug_uID, debug_is_Unitig;
                // get_R_to_U(ruIndex, RrId, &debug_uID, &debug_is_Unitig);
                // if(seed.a[k].t_id != debug_uID) fprintf(stderr, "ERROR1\n");
            }

            q_beg = t_beg = (uint32_t)-1; q_end = t_end = 0;
            for (k = 0; k < seed.n; k++)
            {
                if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;

                if(seed.a[k].q_pos < q_beg) q_beg = seed.a[k].q_pos;
                if(seed.a[k].q_pos > q_end) q_end = seed.a[k].q_pos;
            }

            if(q_beg<=q_end && t_beg<=t_end && tLen > 0 && qLen > 0)
            {
                
                ///exclude extraordinary points of t_pos
                ///note here all q_pos are continual
                ///here we should use base position, instead of the read offset
                if((DIFF((q_end+1-q_beg), (t_end+1-t_beg))) > (q_end+1-q_beg)*extraord_rate)
                {
                    uint32_t leftLen = 0, rightLen = 0, m = 0, median = 0, diff = (q_end+1-q_beg)*(1+extraord_rate)*0.5;
                    for (k = 0; k < seed.n; k++)
                    {
                        median += seed.a[k].t_pos;
                    }
                    median = median / seed.n;

                    leftLen = median - diff;
                    if(median < diff) leftLen = 0;

                    rightLen = median + diff;
                    if(rightLen >= tLen) rightLen = tLen - 1;

                    t_beg = (uint32_t)-1; t_end = 0; m = 0;
                    for (k = 0; k < seed.n; k++)
                    {
                        if(seed.a[k].t_pos >= leftLen && seed.a[k].t_pos <= rightLen)
                        {
                            if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                            if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;
                            seed.a[m] = seed.a[k];
                            m++;
                        }
                    }
                    seed.n = m;

                    if(t_beg>t_end) goto direct_skip;
                }

                ///here we know q_beg, q_end, qLen
                ///and          t_beg, t_end, tLen
                ///there might be two directions: 
                ///a) query and target at the same direction
                ///b) query and target at different direction
                get_contig_overlap_interval(0, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);

                uint32_t flag_forward = get_haplotype_rate(ug, read_g, reverse_sources, ruIndex, b_0, 
                interval_q_beg, interval_q_end, anchor_cId, Hap_rate);
                if(flag_forward)
                {
                    x.t_id = anchor_cId;
                    x.q_pos = interval_q_beg;
                    x.t_pos = interval_q_end;
                    ///x.is_color = 0;
                    x.is_color = tLen;
                    kv_push(Hap_Align, alignment, x);
                }

                get_contig_overlap_interval(1, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);


                uint32_t flag_backward = get_haplotype_rate(ug, read_g, reverse_sources, ruIndex, b_0, 
                interval_q_beg, interval_q_end, anchor_cId, Hap_rate);
                if(flag_backward)
                {
                    x.t_id = anchor_cId;
                    x.q_pos = interval_q_beg;
                    x.t_pos = interval_q_end;
                    ///x.is_color = 1;
                    x.is_color = tLen;
                    kv_push(Hap_Align, alignment, x);
                }
            }
            
        }

        direct_skip:
        i++;
    }

    ///must sort here
    radix_sort_Hap_Align_sort(alignment.a, alignment.a+alignment.n);
    
    uint32_t m = 0;
    ///here is a bug
    sLen = 0; anchor_cId = tLen = q_beg = q_end = (uint32_t)-1;
    for (i = 0; i < alignment.n; i++)
    {
        if(anchor_cId != alignment.a[i].t_id)
        {
            if(anchor_cId != (uint32_t)-1)
            {
                alignment.a[m].t_id = anchor_cId;
                alignment.a[m].q_pos = q_beg;
                alignment.a[m].t_pos = q_end;
                alignment.a[m].is_color = tLen;
                m++;
            }
            anchor_cId = alignment.a[i].t_id;
            sLen = 0;
        }
        if(alignment.a[i].t_pos < alignment.a[i].q_pos) continue;
        if(alignment.a[i].t_pos - alignment.a[i].q_pos + 1 > sLen)
        {
            sLen = alignment.a[i].t_pos - alignment.a[i].q_pos + 1;
            q_beg = alignment.a[i].q_pos;
            q_end = alignment.a[i].t_pos;
            tLen = alignment.a[i].is_color;
        }
    }
    if(anchor_cId != (uint32_t)-1)
    {
        alignment.a[m].t_id = anchor_cId;
        alignment.a[m].q_pos = q_beg;
        alignment.a[m].t_pos = q_end;
        alignment.a[m].is_color = tLen;
        m++;
    }
    alignment.n = m;

    /**********************for debug****************************/
    // fprintf(stderr, "***alignment.m: %u\n", (uint32_t)alignment.n);
    // for (i = 0; i < alignment.n; i++)
    // {
    //     fprintf(stderr, "c_id: %u, beg: %u, end: %u, tLen: %u\n", 
    //      alignment.a[i].t_id, alignment.a[i].q_pos, 
    //     alignment.a[i].t_pos, alignment.a[i].is_color);
    // }
    /**********************for debug****************************/
    asg_arc_t *e;
    for (i = 0; i < alignment.n; i++)
    {
        is_build = 0;
        if(alignment.a[i].t_pos < alignment.a[i].q_pos) continue;
        sLen = alignment.a[i].t_pos - alignment.a[i].q_pos + 1;
        tLen = alignment.a[i].is_color;
        ///if overlap is short, must fully cover one of the read
        if(sLen < long_hap_overlap)
        {
            if(sLen >= (fully_cover_rate*tLen) || sLen >= (fully_cover_rate*qLen))
            {
                is_build = 1;
            }
        }
        else///if is long, can be partly cover one of the read
        {
            if(sLen >= (long_hap_overlap_rate*tLen) || sLen >= (long_hap_overlap_rate*qLen))
            {
                is_build = 1;
            }
        }

        if(is_build)
        {
            e = asg_arc_pushp(bi_g);
            e->del = 0;
            e->ol = sLen;
            e->ul = currentId; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(qLen);
            e->v = alignment.a[i].t_id;

            if(is_bi_edge)
            {
                e = asg_arc_pushp(bi_g);
                e->del = 0;
                e->ol = sLen;
                e->ul = alignment.a[i].t_id; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(qLen);
                e->v = currentId;
            }

        }
    }
    /**
    is_build = 0;
    ///uint32_t long_hap_overlap, float long_hap_overlap_rate
    ///if overlap is short, must fully cover one of the read
    if(sLen < long_hap_overlap)
    {
        if(sLen >= (fully_cover_rate*tLen) || sLen >= (fully_cover_rate*qLen))
        {
            is_build = 1;
        }
    }
    else///if is long, can be partly cover one of the read
    {
        if(sLen >= (long_hap_overlap_rate*tLen) || sLen >= (long_hap_overlap_rate*qLen))
        {
            is_build = 1;
        }
    }

    if(is_build)
    {
        asg_arc_t *e;
        e = asg_arc_pushp(bi_g);
        e->del = 0;
        e->ol = flag;
        e->ul = i; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(u_vecs.a.n);
        e->v = cId;
    }
    **/
    

    /**
    radix_sort_Hap_Align_sort(alignment.a, alignment.a+alignment.n);
    fprintf(stderr, "alignment.n: %u, sLen: %u, q_beg: %u, q_end: %u\n", 
    (uint32_t)alignment.n, sLen, q_beg, q_end);
    for (i = 0; i < alignment.n; i++)
    {
        fprintf(stderr, "dir: %u, c_id: %u, beg: %u, end: %u\n", 
        alignment.a[i].is_color, alignment.a[i].t_id, alignment.a[i].q_pos, 
        alignment.a[i].t_pos);
    }
    **/

    


    /**
    uint32_t self_offset, uId;
    ma_utg_t* reads;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
            {
                fprintf(stderr, "ERROR1\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
                {
                    fprintf(stderr, "ERROR2: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    


    


    iterator.offset = iterator.readI = iterator.untigI = 0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
            {
                fprintf(stderr, "ERROR3\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
                {
                    fprintf(stderr, "ERROR33: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    **/


    kv_destroy(alignment);
    kv_destroy(seed);
    free(b_tid.b.a);
    return 1;
}

inline int get_useful_contig_advance_back(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, kvec_t_u64_warp* u_vecs, buf_t* b_0, uint64_t* contigBeg)
{
    uint32_t i = 0, j, k, sLen, is_found, rId, beg, end;
    uint32_t q_beg, q_end, t_beg, t_end, qLen, tLen;
    uint32_t interval_q_beg, interval_q_end, interval_t_beg, interval_t_end;
    uint64_t anchor_offset, offset;
    uint64_t anchor_cId, cId;

    buf_t b_tid;
    memset(&b_tid, 0, sizeof(buf_t));

    kvec_t(Hap_Align) seed;
    kv_init(seed);
    Hap_Align x;
    rIdContig iterator, iterator_tid;
    iterator.b_0 = b_0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    if(u_vecs->a.n == 0) return -1;

    qLen = get_contig_len(ug, b_0);
    
    /**********************for debug****************************/
    // if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
    // {
    //     i = 0;
    //     fprintf(stderr,"*******\n");
    //     for (i = 0; i < u_vecs->a.n; i++)
    //     {
    //         offset = u_vecs->a.a[i]>>33; 
    //         cId = (uint32_t)(u_vecs->a.a[i]);
    //         fprintf(stderr, "((%u) off: %u, cId: %d, HAP: %u)\n", i, offset, (int)cId, 
    //         ((u_vecs->a.a[i]&(uint64_t)(0x100000000))!=0));
    //     }
    // }
    /**********************for debug****************************/


    i = 0;
    while (i < u_vecs->a.n)
    {
        ///there are three cases:
        /**
         * 1) u_vecs->a.a[i] has haplotype infor at the primary contigs
         * 2) u_vecs->a.a[i] has haplotype infor at the alternative contigs, 
         *    (uint32_t)u_vecs->a.a[i]==(uint32_t)-1
         * 3) u_vecs->a.a[i] itself is labled with HAP_LABLE, 
         *    u_vecs->a.a[i]&(uint64_t)(0x100000000)>0 && (uint32_t)u_vecs->a.a[i]==(uint32_t)-1
         * 4) u_vecs->a.a[i] does not have any haplotype infor. In this case, 
         *    self_offsetLen is not consecutive 
         * 
        **/
        if((u_vecs->a.a[i]&(uint64_t)(0x100000000)) || ((uint32_t)u_vecs->a.a[i]==((uint32_t)-1)))
        {
            i++;
            continue;
        }
        

        anchor_offset = (u_vecs->a.a[i]>>33); anchor_cId = (uint32_t)(u_vecs->a.a[i]); 
        sLen = 0; is_found = 0; seed.n = 0;
        for (j = i; j < u_vecs->a.n; j++)
        {
            ///offset is the self offset at query
            ///cId is the target Id (contig)
            offset = u_vecs->a.a[j]>>33; 
            cId = (uint32_t)(u_vecs->a.a[j]);

            ///means we go ahead one step at query
            if(offset != anchor_offset)
            {
                ///sLen must be >=1
                if(is_found == 0) 
                {
                    break;
                }
                is_found = 0;
                sLen++;
                anchor_offset = offset;
            } 

            if(cId == anchor_cId)
            {
                x.q_pos = u_vecs->a.a[j]>>33;
                x.t_id = (uint32_t)(u_vecs->a.a[j]);
                x.is_color = 1;
                x.t_pos = (uint32_t)-1;
                u_vecs->a.a[j] = u_vecs->a.a[j]|0xffffffff;
                is_found++;
                // if(seed.n > 0 && seed.a[seed.n-1].t_id == x.t_id && 
                // seed.a[seed.n-1].q_pos == x.q_pos)
                // {
                //     seed.a[seed.n-1].is_color++;
                //     continue;
                // }
                kv_push(Hap_Align, seed, x);
            } 
        }

        if(j == u_vecs->a.n && is_found > 0)
        {
            sLen++;
        }


        // for (j = 1; j < seed.n; j++)
        // {
        //     if(seed.a[j].t_id != seed.a[j-1].t_id)
        //     {
        //         fprintf(stderr, "hehe\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }

        //     if(seed.a[j].q_pos < seed.a[j-1].q_pos)
        //     {
        //         fprintf(stderr, "haha\n");
        //         fprintf(stderr, "sLen: %u, seed.n: %u\n", sLen, seed.n);
        //     }
        // }

        ///don't need to know how long of this seed
        if(sLen >= hap_seed && seed.n > 0)
        {
            /**********************for debug****************************/
            ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
            ///fprintf(stderr, "\nsLen: %u, seed.n: %u, anchor_cId: %u\n", sLen, seed.n, anchor_cId);
            /**********************for debug****************************/
            
            uint32_t inner_off = 0, pre_q_pos = (uint32_t)-1, RrId;
            iterator.offset = iterator.readI = iterator.untigI = 0;

            ///for one vector, all t_id should be same; that is why we recover tid here
            beg = (uint32_t)(contigBeg[seed.a[0].t_id]);
            b_tid.b.n = 0;ug->u.a[beg>>1].circ = 0;
            get_long_tip_length(ug->g, &(ug->u), beg, &end, &b_tid);
            iterator_tid.b_0 = &b_tid;
            iterator_tid.offset = iterator_tid.readI = iterator_tid.untigI = 0;
            ///we have already got qLen at the begining
            tLen = get_contig_len(ug, &b_tid);

            for (k = 0; k < seed.n; k++)
            {

                rId = get_rId_from_contig_by_offset(ug, &iterator, seed.a[k].q_pos);
                ///actually we shouldn't have this case
                if(rId == (uint32_t)-1) continue;
                if(pre_q_pos != seed.a[k].q_pos)
                {
                    inner_off = 0;
                    pre_q_pos = seed.a[k].q_pos;
                }
                else
                {
                    inner_off++;
                }
                
                RrId = get_reverseId(read_g, reverse_sources, ruIndex, rId, seed.a[k].t_id, 
                inner_off);
                ///if(RrId == (uint32_t)-1) fprintf(stderr, "ERROR1\n");
                ///for one vector, all t_id should be same
                seed.a[k].t_pos = get_offset_from_contig_by_rId(ug, &iterator_tid, RrId);
                // if(seed.a[k].t_pos == (uint32_t)-1) fprintf(stderr, "ERROR2\n");
                // uint32_t debug_uID, debug_is_Unitig;
                // get_R_to_U(ruIndex, RrId, &debug_uID, &debug_is_Unitig);
                // if(seed.a[k].t_id != debug_uID) fprintf(stderr, "ERROR1\n");
                
                /**********************for debug****************************/
                ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
                // fprintf(stderr, "%u, q_pos: %u, t_id: %u, t_pos: %u, inner_off: %u, rId: %u, RrId: %u\n", 
                // k, seed.a[k].q_pos, seed.a[k].t_id, seed.a[k].t_pos, inner_off, rId, RrId);
                /**********************for debug****************************/
            }

            q_beg = t_beg = (uint32_t)-1; q_end = t_end = 0;
            for (k = 0; k < seed.n; k++)
            {
                if(seed.a[k].t_pos < t_beg) t_beg = seed.a[k].t_pos;
                if(seed.a[k].t_pos > t_end) t_end = seed.a[k].t_pos;

                if(seed.a[k].q_pos < q_beg) q_beg = seed.a[k].q_pos;
                if(seed.a[k].q_pos > q_end) q_end = seed.a[k].q_pos;
            }

            

            if(q_beg<=q_end && t_beg<=t_end)
            {
                ///here we know q_beg, q_end, qLen
                ///and          t_beg, t_end, tLen
                ///there might be two directions: 
                ///a) query and target at the same direction
                ///b) query and target at different direction
                get_contig_overlap_interval(0, q_beg, q_end, qLen, 
                t_beg, t_end, tLen, &interval_q_beg, &interval_q_end,
                &interval_t_beg, &interval_t_end);


                /**********************for debug****************************/
                ///if(u_vecs->a.n > 20 && u_vecs->a.n < 50)
                {
                    fprintf(stderr, "q_beg: %u, q_end: %u, qLen: %u, t_beg: %u, t_end: %u, tLen: %u\n", 
                    q_beg, q_end, qLen, t_beg, t_end, tLen);

                    fprintf(stderr, "interval_q_beg: %u, interval_q_end: %u, interval_t_beg: %u, interval_t_end: %u\n", 
                    interval_q_beg, interval_q_end, interval_t_beg, interval_t_end);   
                }
                /**********************for debug****************************/
            
            }
            
        }
        i++;
    }






    /**
    uint32_t self_offset, uId;
    ma_utg_t* reads;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
            {
                fprintf(stderr, "ERROR1\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_rId_from_contig_by_offset(ug, &iterator, self_offset)!=rId)
                {
                    fprintf(stderr, "ERROR2: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    


    


    iterator.offset = iterator.readI = iterator.untigI = 0;
    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (j = 0, self_offset = 0; j < b_0->b.n; j++)
    {
        uId = b_0->b.a[j]>>1;
        ///if(IsMerge(ug->u, uId)>0) continue;
        reads = &(ug->u.a[uId]);
        ///scan all reads
        ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
        for (k = 0; k < reads->n; k++, self_offset++)
        {
            rId = reads->a[k]>>33;
            if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
            {
                fprintf(stderr, "ERROR3\n");
            }
        }
    }

    iterator.offset = iterator.readI = iterator.untigI = 0;
    for (long long debug_i = self_offset; debug_i >= 0; debug_i--)
    {
        for (j = 0, self_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                if(debug_i != self_offset) continue;
                rId = reads->a[k]>>33;
                if(get_offset_from_contig_by_rId(ug, &iterator, rId)!=self_offset)
                {
                    fprintf(stderr, "ERROR33: self_offset: %u\n", self_offset);
                }
            }
        }
    }
    **/



    kv_destroy(seed);
    free(b_tid.b.a);
    return 1;
}


#define contig_seed 20
inline int get_useful_contig(kvec_t_u64_warp* u_vecs, float density, uint32_t miniLen, uint32_t* r_cId)
{
    (*r_cId) = (uint32_t)-1;
    uint32_t cId;
    uint32_t intervalLen = 0, realLen = 0, useful_index = (uint32_t)-1;
    uint32_t i = u_vecs->i, end;
    float T_density = density/2;
    
    if(i >= u_vecs->a.n) return -1;
    if((uint32_t)(u_vecs->a.a[i]) == (uint32_t)-1)
    {
        u_vecs->i++;
        return 0;
    } 

    end = i + contig_seed;
    if(end > u_vecs->a.n) end = u_vecs->a.n;
    cId = (uint32_t)(u_vecs->a.a[i]);
    for (; i < end; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i])) realLen++;
        intervalLen++;
        if(realLen >= intervalLen*density) useful_index = i;
    }

    if(realLen < intervalLen*T_density)
    {
        u_vecs->i++;
        return 0;
    } 
    ///here we found a useful seed
    ///it seems we don't need to scan backward
    for (; i < u_vecs->a.n; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i])) realLen++;
        intervalLen++;
        if(realLen < intervalLen*T_density) break;
        if(realLen >= intervalLen*density) useful_index = i;
    }

    if(useful_index == (uint32_t)-1 || (useful_index - u_vecs->i) < miniLen)
    {
        u_vecs->i++;
        return 0;
    }

    ///don't need to set backward
    end = (uint32_t)-1;
    for (i = u_vecs->i; i < u_vecs->a.n; i++)
    {
        if(cId == (uint32_t)(u_vecs->a.a[i]))
        {
            ///u_vecs->a.a[i] = (uint32_t)-1;
            u_vecs->a.a[i] = u_vecs->a.a[i]|(uint64_t)(0xffffffff);
        } 
    }
    (*r_cId) = cId;
    u_vecs->i++;
    return (useful_index - u_vecs->i);
}





#define UNVISIT (uint32_t)(0x7fffffff)
#define RED 0
#define BLACK 1
#define ISO 2
#define LABLE 3
void bi_paration(asg_t *bi_g, uint64_t* array, uint32_t bi_graph_Len)
{
    
    kvec_t(uint64_t) a;
    kv_init(a);
    a.a = array;
    a.n = a.m = bi_g->n_seq;
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);
    uint32_t i, len, v, w, k, nv, roundID = 0;
    asg_arc_t *av;
    for (i = 0; i < bi_g->n_seq; i++)
    {
        bi_g->seq[i].c = 0;
    }

    re_partition:

    for (i = 0; i < bi_g->n_seq; i++)
    {
        /****************************may have bugs********************************/
        ///how many reads contained in this contig
        len = (a.a[i]>>33);
        /****************************may have bugs********************************/
        ///if the contig is too small, or the contig has already been visited
        ///if(len < bi_graph_Len || bi_g->seq[i].len != UNVISIT) continue;
        if(len < bi_graph_Len || bi_g->seq[i].c == 1) continue;
        if(roundID == 0 && bi_g->seq[i].len == UNVISIT) continue;

        ///set the color of this node
        bi_g->seq[i].len = RED; bi_g->seq[i].c = 1;
        kdq_push(uint32_t, buf, i);
        while (kdq_size(buf) != 0)
        {
            v = *(kdq_pop(uint32_t, buf)); bi_g->seq[v].c = 1;
            if(bi_g->seq[v].len == ISO) continue;
            nv = asg_arc_n(bi_g, v);
            av = asg_arc_a(bi_g, v);


            ///get all out-nodes of v
            for(k = 0; k < nv; k++)
            {
                w = av[k].v;
                /****************************may have bugs********************************/
                ///if(bi_g->seq[w].len == bi_g->seq[v].len)
                len = (a.a[w]>>33);
                ///only check large contig
                if(len >= bi_graph_Len && bi_g->seq[w].len == bi_g->seq[v].len)
                { /****************************may have bugs********************************/
                    break;
                }
            }
            ///means v is conflict 
            if(k != nv)
            {
                bi_g->seq[v].len = ISO; bi_g->seq[v].c = 1;
                continue;
            }



            ///if v is not conflict with others
            for(k = 0; k < nv; k++)
            {
                w = av[k].v;
                ///if the out-node has not been visited
                if(bi_g->seq[w].len == UNVISIT)
                {
                    bi_g->seq[w].len = 1 - bi_g->seq[v].len;
                    bi_g->seq[w].c = 1;
                    /****************************may have bugs********************************/
                    len = (a.a[w]>>33);
                    /****************************may have bugs********************************/
                    if(len < bi_graph_Len) continue;
                    kdq_push(uint32_t, buf, w);
                }
                ///here bi_g->seq[w].len might be ISO or another color
                ///don't need to do anything here
            }
        }
    }

    if(roundID == 0)
    {
        roundID = 1;
        goto re_partition;
    }
    

    //secondary checking
    // uint32_t a_color;
    // kdq_size(buf) = 0;
    // for (i = 0; i < bi_g->n_seq; i++)
    // {
    //     /****************************may have bugs********************************/
    //     len = (a.a[i]>>33);
    //     /****************************may have bugs********************************/
    //     ///just check end contig
    //     if(bi_g->seq[i].len != UNVISIT || (a.a[i] & (uint64_t)(0x100000000)) == 0) continue;
        
    //     v = i;
    //     nv = asg_arc_n(bi_g, v);
    //     av = asg_arc_a(bi_g, v);
    //     a_color = UNVISIT;

    //     for(k = 0; k < nv; k++)
    //     {
    //         w = av[k].v;
    //         ///check all out-nodes that has already been colored
    //         if(bi_g->seq[w].len != UNVISIT)
    //         {
    //             ///if this is the first colored node
    //             if(a_color == UNVISIT)
    //             {
    //                 a_color = bi_g->seq[w].len;
    //             }///if this is not
    //             else if(a_color != bi_g->seq[w].len)
    //             {
    //                 a_color = ISO;
    //                 bi_g->seq[v].len = ISO;
    //                 break;
    //             }
    //         }
    //     }

    //     if(a_color == ISO) continue;
    //     ///no colored out-node
    //     if(a_color == UNVISIT)
    //     {
    //         bi_g->seq[v].len = RED;
    //     }
    //     else
    //     {
    //         bi_g->seq[v].len = 1 - a_color;
    //     }



        

    //     ///check if all reachable nodes are end-contig
    //     kdq_push(uint32_t, buf, i);
    //     while (kdq_size(buf) != 0)
    //     {
    //         v = *(kdq_pop(uint32_t, buf));
    //         ///find a non-end contig
    //         if((a.a[v] & (uint64_t)(0x100000000)) == 0)
    //         {
    //             bi_g->seq[i].len = UNVISIT;
    //             kdq_size(buf) = 0;
    //             break;
    //         }

    //         nv = asg_arc_n(bi_g, v);
    //         av = asg_arc_a(bi_g, v);
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             ///if the out-node has not been visited
    //             if(bi_g->seq[w].len == UNVISIT)
    //             {
    //                 bi_g->seq[v].len = LABLE;
    //                 kdq_push(uint32_t, buf, w);
    //             }
    //         }
    //     }

    //     for (v = 0; v < bi_g->n_seq; v++)
    //     {
    //         if(bi_g->seq[v].len == LABLE) bi_g->seq[v].len = UNVISIT;
    //     }

    //     if(bi_g->seq[i].len == UNVISIT)
    //     {
    //         continue;
    //     }






    //     ///now all reachable nodes of i are end-contigs
    //     kdq_push(uint32_t, buf, i);
    //     while (kdq_size(buf) != 0)
    //     {
    //         ///each v here is the uncolored node in the first round
    //         v = *(kdq_pop(uint32_t, buf));
    //         if(bi_g->seq[v].len == ISO) continue;
    //         nv = asg_arc_n(bi_g, v);
    //         av = asg_arc_a(bi_g, v);


    //         ///get all out-nodes of v
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             if(bi_g->seq[w].len == bi_g->seq[v].len)
    //             {
    //                 break;
    //             }
    //         }
    //         ///means v is conflict 
    //         if(k != nv)
    //         {
    //             bi_g->seq[v].len = ISO;
    //             continue;
    //         }



    //         ///if v is not conflict with others
    //         for(k = 0; k < nv; k++)
    //         {
    //             w = av[k].v;
    //             ///if the out-node has not been visited
    //             if(bi_g->seq[w].len == UNVISIT)
    //             {   
    //                 bi_g->seq[w].len = 1 - bi_g->seq[v].len;
    //                 kdq_push(uint32_t, buf, w);
    //             }
    //         }
    //     }
    // }

    kdq_destroy(uint32_t, buf);
}

void further_clean_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniLen, uint32_t bi_graph_Len,
uint32_t long_hap_overlap, float long_hap_overlap_rate, float lable_match_rate)
{
    asg_t *bi_g = NULL;
    bi_g = asg_init();
    kvec_t(uint64_t) a;
    kv_init(a);
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, beg = 0, end, uId, cId, rId, i, j, k, self_offset, self_label_offset;
    uint64_t uInfor = 0;
    ma_utg_t* reads;
    ///int flag;
    memset(visit, 0, nsg->n_seq);
    /******************************set ruIndex*********************************/
    for (v = 0; v < n_vtx; ++v) 
    {
        beg = v;
        if(nsg->seq[v>>1].c == ALTER_LABLE || nsg->seq[beg>>1].del || visit[beg>>1])
        {
            continue;
        }

        ///why I have this line? It might be wrong
        ///for example, if num(beg) == 0, and num(beg^1) == 2, beg is a unitig
        if(get_real_length(nsg, beg, NULL)<=0 && get_real_length(nsg, beg^1, NULL)>0)
        {
            continue;
        }


        if(get_real_length(nsg, beg^1, NULL) == 1)
        {
            get_real_length(nsg, beg^1, &end);
            if(get_real_length(nsg, end^1, NULL) == 1)
            {
                continue;
            }
        }

        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);

        uInfor = 0;
        //scan all untigs
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            visit[uId] = 1;
            reads = &(ug->u.a[uId]);
            uInfor += reads->n;
        }
        /****************************may have bugs********************************/
        ///uInfor = uInfor << 32; uInfor = uInfor | (uint64_t)beg;
        uInfor = uInfor << 33; uInfor = uInfor | (uint64_t)beg;
        /****************************may have bugs********************************/
        kv_push(uint64_t, a, uInfor);
    }

    ///sort by number of reads in a contig
    radix_sort_arch64(a.a, a.a + a.n);
    for (i = 0; i < (a.n>>1); ++i) 
    {
        uInfor = a.a[i];
        a.a[i] = a.a[a.n - i - 1];
        a.a[a.n - i - 1] = uInfor;
    }


    for (v = 0; v < a.n; v++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[v]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        cId = v;
        ///separated contig
        if(get_real_length(nsg, beg^1, NULL) == 0 && get_real_length(nsg, end, NULL) == 0)
        {
            a.a[v] = a.a[v] | (uint64_t)(0x100000000);
        }
        ///set the contig Id for each read
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            reads = &(ug->u.a[uId]);
            for (j = 0; j < reads->n; j++)
            {
                rId = reads->a[j]>>33;
                set_R_to_U(ruIndex, rId, cId, 1);
            }
        }
    }
    /******************************set ruIndex*********************************/

    
    ///scan all contigs from longest one to the shortest one
    for (i = 0; i < a.n; i++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        
        
        u_vecs.a.n = 0;
        ///scan all untigs of this contig
        for (j = 0, self_offset = 0, self_label_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                rId = reads->a[k]>>33;
                ///if(read_g->seq[rId].c == HAP_LABLE) continue;
                ///self_offset is the offset of this read in contig
                query_reverse_sources(read_g, reverse_sources, ruIndex, rId, self_offset, &u_vecs);

                if(read_g->seq[rId].c == HAP_LABLE) self_label_offset++;
            }
        }

        if(self_offset >= long_hap_overlap && self_label_offset >= (self_offset*lable_match_rate))
        {
            asg_seq_set(bi_g, i, RED, 0);
        }
        else
        {
            asg_seq_set(bi_g, i, UNVISIT, 0);
        }

        ///fprintf(stderr, "self_label_offset: %u, self_offset: %u\n", self_label_offset, self_offset);
        
        get_useful_contig_advance(ug, read_g, reverse_sources, ruIndex, &u_vecs, b_0, a.a, 
        bi_g, i, density, long_hap_overlap, long_hap_overlap_rate, 1);
    }

    asg_cleanup(bi_g);

    ///fprintf(stderr, "***********n_seq: %u, n_arc: %u\n", bi_g->n_seq, bi_g->n_arc);
    for (v = 0; v < bi_g->n_seq; v++)
    {
        uint32_t nv = asg_arc_n(bi_g, v), nw, w;
        asg_arc_t *av = asg_arc_a(bi_g, v), *aw;
        for (i = 0; i < nv; ++i)
        {
            w = av[i].v;
            if(w == v)
            {
                av[i].del = 1;
                continue;
            }
            nw = asg_arc_n(bi_g, w);
            aw = asg_arc_a(bi_g, w);
            for (j = 0; j < nw; j++)
            {
                if(aw[j].v == v) break;
            }

            if(j == nw) av[i].del = 1;
        }
    }

    asg_cleanup(bi_g);
    ///fprintf(stderr, "***********n_seq: %u, n_arc: %u\n", bi_g->n_seq, bi_g->n_arc);

    bi_paration(bi_g, a.a, bi_graph_Len);

    for (v = 0; v < bi_g->n_seq; v++)
    {
        beg = (uint32_t)(a.a[v]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);

        if(bi_g->seq[v].len == BLACK)
        {
            for (i = 0; i < b_0->b.n; i++)
            {
                nsg->seq[(b_0->b.a[i])>>1].c = ALTER_LABLE;
            }

            for (i = 0; i < b_0->b.n; i++)
            {
                asg_seq_drop(nsg, ((b_0->b.a[i])>>1));
            }
        }
        







        // fprintf(stderr, "cId: %u, beg>>1: %u, end>>1: %u, b_0->b.n: %u, Len: %u, type: %u\n", 
        // v, beg>>1, end>>1, (uint32_t)b_0->b.n, (uint32_t)(a.a[v]>>33), bi_g->seq[v].len);
        // uint32_t nv = asg_arc_n(bi_g, v), w;
        // asg_arc_t *av = asg_arc_a(bi_g, v);
        // for (i = 0; i < nv; ++i)
        // {
        //     w = av[i].v;
        //     fprintf(stderr, "w: %u\n", w);
        // }
    }

    
    /**
    for (i = 1; i < a.n; i++)
    {
        if((a.a[i]>>33) > (a.a[i-1]>>33))
        {
            fprintf(stderr, "hehe\n");
        }
    }

    
    for (i = 0; i < a.n; i++)
    {
        uint32_t k, get_cId, tLen = 0, is_Unitig;
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        for (j = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            reads = &(ug->u.a[uId]);
            tLen += reads->n;
            for (k = 0; k < reads->n; k++)
            {
                rId = reads->a[k]>>33;
                get_R_to_U(ruIndex, rId, &get_cId, &is_Unitig);
                if(is_Unitig != 1 || get_cId != i)
                {
                    fprintf(stderr, "###is_Unitig: %u, get_cId: %u, i: %u\n", 
                    is_Unitig, get_cId, i);
                }
            }
        }

        uint32_t qLen = 0, m;
        for (m = 0; m < ruIndex->len; m++)
        {
            get_R_to_U(ruIndex, m, &get_cId, &is_Unitig);
            if(get_cId == (uint32_t)-1 || is_Unitig != 1) continue;
            if(get_cId == i)
            {
                qLen = 0;
                for (j = 0; j < b_0->b.n; j++)
                {
                    uId = b_0->b.a[j]>>1;
                    reads = &(ug->u.a[uId]);
                    tLen += reads->n;
                    for (k = 0; k < reads->n; k++)
                    {
                        rId = reads->a[k]>>33;
                        if(m == rId)
                        {
                            qLen = 1;
                            goto found;
                        }
                    }
                }
                found:;
                if(qLen == 0)
                {
                    fprintf(stderr, "***is_Unitig: %u, get_cId: %u, m: %u\n", 
                    is_Unitig, get_cId, m);
                }
            }
        }
    }
    **/

    asg_cleanup(nsg);
    free(a.a);
    kv_destroy(u_vecs.a);
    asg_destroy(bi_g);
}


void process_bi_graph(asg_t *bi_g)
{
    asg_cleanup(bi_g);
    asg_arc_del_multi(bi_g);
	bi_g->is_symm = 1;
    uint32_t v, i, j;
    for (v = 0; v < bi_g->n_seq; v++)
    {
        uint32_t nv = asg_arc_n(bi_g, v), nw, w;
        asg_arc_t *av = asg_arc_a(bi_g, v), *aw;
        for (i = 0; i < nv; ++i)
        {
            w = av[i].v;
            if(w == v)
            {
                av[i].del = 1;
                continue;
            }
            nw = asg_arc_n(bi_g, w);
            aw = asg_arc_a(bi_g, w);
            for (j = 0; j < nw; j++)
            {
                if(aw[j].v == v) break;
            }

            if(j == nw) av[i].del = 1;
        }
    }

    asg_cleanup(bi_g);

    ///asg_symm(bi_g);
}

void further_clean_untig_graph_trio(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t bi_graph_Len,
uint32_t long_hap_overlap, float long_hap_overlap_rate, float lable_match_rate)
{
    asg_t *bi_g = NULL;
    bi_g = asg_init();
    kvec_t(uint64_t) a;
    kv_init(a);
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, beg = 0, end, uId, cId, rId, i, j, k, self_offset, self_label_offset;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    uint64_t uInfor = 0;
    ma_utg_t* reads;
    ///int flag;
    memset(visit, 0, nsg->n_seq);


    /****************************may have bugs********************************/
    for (v = 0; v < nsg->n_seq; ++v) 
    {
        uId = v;
        if(nsg->seq[uId].c != HAP_LABLE) continue;
        reads = &(ug->u.a[uId]);
        for (k = 0; k < reads->n; k++)
        {
            rId = reads->a[k]>>33;
            read_g->seq[rId].c = HAP_LABLE;
        }
    }
    /****************************may have bugs********************************/



    /******************************set ruIndex*********************************/
    for (v = 0; v < n_vtx; ++v) 
    {
        beg = v;
        if(nsg->seq[v>>1].c == ALTER_LABLE || nsg->seq[beg>>1].del || visit[beg>>1])
        {
            continue;
        }

        if(get_real_length(nsg, beg^1, NULL) == 1)
        {
            get_real_length(nsg, beg^1, &end);
            if(get_real_length(nsg, end^1, NULL) == 1) continue;
        }

        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);

        uInfor = 0;
        //scan all untigs
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            visit[uId] = 1;
            reads = &(ug->u.a[uId]);
            uInfor += reads->n;
        }
        /****************************may have bugs********************************/
        ///uInfor = uInfor << 32; uInfor = uInfor | (uint64_t)beg;
        uInfor = uInfor << 33; uInfor = uInfor | (uint64_t)beg;
        /****************************may have bugs********************************/
        kv_push(uint64_t, a, uInfor);
    }

    ///sort by number of reads in a contig
    radix_sort_arch64(a.a, a.a + a.n);
    for (i = 0; i < (a.n>>1); ++i) 
    {
        uInfor = a.a[i];
        a.a[i] = a.a[a.n - i - 1];
        a.a[a.n - i - 1] = uInfor;
    }


    for (v = 0; v < a.n; v++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[v]);

        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);
       
        cId = v;
        ///individual contig
        if(get_real_length(nsg, beg^1, NULL) == 0 && get_real_length(nsg, end, NULL) == 0)
        {
            a.a[v] = a.a[v] | (uint64_t)(0x100000000);
        }
        ///set the contig Id for each read
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            reads = &(ug->u.a[uId]);
            for (j = 0; j < reads->n; j++)
            {
                rId = reads->a[j]>>33;
                set_R_to_U(ruIndex, rId, cId, 1);
            }
        }
    }
    /******************************set ruIndex*********************************/

    
    ///scan all contigs from the longest one to the shortest one
    for (i = 0; i < a.n; i++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);

        u_vecs.a.n = 0;
        ///scan all unitigs of this contig
        for (j = 0, self_offset = 0, self_label_offset = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            ///if(IsMerge(ug->u, uId)>0) continue;
            reads = &(ug->u.a[uId]);
            ///scan all reads
            ///self_offset will skip fake(merge) nodes, but not skip HAP_LABLE
            for (k = 0; k < reads->n; k++, self_offset++)
            {
                rId = reads->a[k]>>33;
                ///if(read_g->seq[rId].c == HAP_LABLE) continue;
                ///self_offset is the offset of this read in contig
                query_reverse_sources(read_g, reverse_sources, ruIndex, rId, self_offset, &u_vecs);

                if(read_g->seq[rId].c == HAP_LABLE) self_label_offset++;
            }
        }

        if(self_offset >= long_hap_overlap && self_label_offset >= (self_offset*lable_match_rate))
        {
            asg_seq_set(bi_g, i, RED, 0);
        }
        else
        {
            asg_seq_set(bi_g, i, UNVISIT, 0);
        }


        ///a.a saves all offest and its corresponding contig ID
        get_useful_contig_advance(ug, read_g, reverse_sources, ruIndex, &u_vecs, b_0, a.a, 
        bi_g, i, density, long_hap_overlap, long_hap_overlap_rate, 1);
    }

    process_bi_graph(bi_g);

    bi_paration(bi_g, a.a, bi_graph_Len);

    for (v = 0; v < bi_g->n_seq; v++)
    {
        beg = (uint32_t)(a.a[v]);
        /**
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        **/
        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);

        if(bi_g->seq[v].len == BLACK)
        {
            for (i = 0; i < b_0->b.n; i++)
            {
                nsg->seq[(b_0->b.a[i])>>1].c = ALTER_LABLE;
            }

            for (i = 0; i < b_0->b.n; i++)
            {
                asg_seq_drop(nsg, ((b_0->b.a[i])>>1));
            }
        }
        







        // fprintf(stderr, "cId: %u, beg>>1: %u, end>>1: %u, b_0->b.n: %u, Len: %u, type: %u\n", 
        // v, beg>>1, end>>1, (uint32_t)b_0->b.n, (uint32_t)(a.a[v]>>33), bi_g->seq[v].len);
        // uint32_t nv = asg_arc_n(bi_g, v), w;
        // asg_arc_t *av = asg_arc_a(bi_g, v);
        // for (i = 0; i < nv; ++i)
        // {
        //     w = av[i].v;
        //     fprintf(stderr, "w: %u\n", w);
        // }
    }

    
    /**
    for (i = 1; i < a.n; i++)
    {
        if((a.a[i]>>33) > (a.a[i-1]>>33))
        {
            fprintf(stderr, "hehe\n");
        }
    }

    
    for (i = 0; i < a.n; i++)
    {
        uint32_t k, get_cId, tLen = 0, is_Unitig;
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        for (j = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            reads = &(ug->u.a[uId]);
            tLen += reads->n;
            for (k = 0; k < reads->n; k++)
            {
                rId = reads->a[k]>>33;
                get_R_to_U(ruIndex, rId, &get_cId, &is_Unitig);
                if(is_Unitig != 1 || get_cId != i)
                {
                    fprintf(stderr, "###is_Unitig: %u, get_cId: %u, i: %u\n", 
                    is_Unitig, get_cId, i);
                }
            }
        }

        uint32_t qLen = 0, m;
        for (m = 0; m < ruIndex->len; m++)
        {
            get_R_to_U(ruIndex, m, &get_cId, &is_Unitig);
            if(get_cId == (uint32_t)-1 || is_Unitig != 1) continue;
            if(get_cId == i)
            {
                qLen = 0;
                for (j = 0; j < b_0->b.n; j++)
                {
                    uId = b_0->b.a[j]>>1;
                    reads = &(ug->u.a[uId]);
                    tLen += reads->n;
                    for (k = 0; k < reads->n; k++)
                    {
                        rId = reads->a[k]>>33;
                        if(m == rId)
                        {
                            qLen = 1;
                            goto found;
                        }
                    }
                }
                found:;
                if(qLen == 0)
                {
                    fprintf(stderr, "***is_Unitig: %u, get_cId: %u, m: %u\n", 
                    is_Unitig, get_cId, m);
                }
            }
        }
    }
    **/

   uint32_t is_Unitig;
   for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    asg_cleanup(nsg);
    free(a.a);
    kv_destroy(u_vecs.a);
    asg_destroy(bi_g);
}

inline void reset_visit_flag(uint8_t* visit, asg_t *read_g, R_to_U* ruIndex, uint32_t contigNum, 
ma_hit_t_alloc* x)
{
    uint32_t k, rId, is_Unitig, Hap_cId;

    if(x->length*2 > contigNum)
    {
        memset(visit, 0, contigNum);
    }
    else
    {
        for (k = 0; k < x->length; k++)
        {
            rId = Get_tn(x->buffer[k]);
        
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }

            ///there are two cases: 
            ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
            ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
            get_R_to_U(ruIndex, rId, &Hap_cId, &is_Unitig);
            if(is_Unitig == 0 || Hap_cId == (uint32_t)-1) continue;
            ///here rId is the id of the read coming from the different haplotype
            ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
            visit[Hap_cId] = 0;
        }
    }
}

void debug_visit_flag(uint8_t* visit, uint32_t contigNum)
{
    uint32_t k;
    for (k = 0; k < contigNum; k++)
    {
        if(visit[k] != 0) fprintf(stderr, "[E::%s]ERROR visit\n", __func__);
    }
}


uint32_t get_readSeq(ma_ug_t *ug, asg_t *read_g, kvec_t_u32_warp* x_vecs, uint64_t* cBeg, 
buf_t* b_0, uint32_t cId)
{
    ma_utg_t* reads = NULL;
    uint32_t beg, end, i, j, rId, xLen, uId, uOri;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    uint64_t tmp;

    x_vecs->a.n = 0;
    beg = (uint32_t)(cBeg[cId]);
    b_0->b.n = 0;
    get_unitig(ug->g, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);
    xLen = cBeg[cId]>>33;
    kv_resize(uint32_t, x_vecs->a, xLen);

    for (i = 0; i < b_0->b.n; i++)
    {
        uId = b_0->b.a[i]>>1;
        uOri = b_0->b.a[i]&(uint32_t)1;
        reads = &(ug->u.a[uId]);
        for (j = 0; j < reads->n; j++)
        {
            if(uOri == 1)
            {
                rId = reads->a[reads->n - j - 1]>>33;
            }
            else
            {
                rId = reads->a[j]>>33;
            }


            tmp = rId<<1;
            if(read_g->seq[rId].c == HAP_LABLE)
            {
                tmp = tmp | 1;
                kv_push(uint32_t, x_vecs->a, tmp);
                continue;
            }
            kv_push(uint32_t, x_vecs->a, tmp);
        }
    }



    if(x_vecs->a.n != xLen) fprintf(stderr, "[E::%s]ERROR: different length\n", __func__);
    return xLen;
}


uint32_t inline retrieve_skip_overlaps(ma_hit_t_alloc* x, uint32_t target)
{
    if(x == NULL) return (uint32_t)-1;

    uint32_t i;
    for (i = 0; i < x->length; i++)
    {
        if(Get_tn(x->buffer[i]) == target)
        {
            return i;
        }
    }

    return (uint32_t)-1;
}

inline uint32_t check_duplicate(Hap_Align_warp* u_buffer, uint32_t x_pos, uint32_t y_pos)
{
    if(u_buffer->x.n == 0) return 0;

    int i = u_buffer->x.n;
    for (i--; i >= 0; i--)
    {
        if(x_pos != (uint32_t)-1 && u_buffer->x.a[i].q_pos != x_pos) return 0;
        if(y_pos != (uint32_t)-1 && u_buffer->x.a[i].t_pos == y_pos)
        {
            u_buffer->x.a[i].is_color++;
            return 1;
        } 
    }

    if(x_pos == (uint32_t)-1) fprintf(stderr, "[E::%s]ERROR: cannot found y\n", __func__);
    
    return 0;
}

///x_vecs->a.a[k]>>1

void get_hap_similarity(uint32_t* list, uint32_t Len, uint32_t target_uId, 
ma_hit_t_alloc* reverse_sources, asg_t *read_g, R_to_U* ruIndex, 
double* Match, double* Total)
{
    #define CUTOFF_THRES 100
    uint32_t i, j, qn, tn, is_Unitig, uId, min_count = 0, max_count = 0, cutoff = 0;;
    for (i = 0; i < Len; i++)
    {
        if(cutoff > CUTOFF_THRES)
        {
            max_count = 0;
            min_count = Len;
        }
        qn = list[i]>>1;
        if(reverse_sources[qn].length > 0) min_count++;
        for (j = 0; j < reverse_sources[qn].length; j++)
        {
            tn = Get_tn(reverse_sources[qn].buffer[j]);
            if(read_g->seq[tn].del == 1)
            {
                get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
            }


            get_R_to_U(ruIndex, tn, &uId, &is_Unitig);
            if(uId!=(uint32_t)-1 && is_Unitig == 1 && uId == target_uId)
            {
                max_count++;
                break;
            }
        }

        //means no match
        if(j == reverse_sources[qn].length)
        {
            cutoff++;
        }
        else
        {
            cutoff = 0;
        }
    }

    (*Match) = max_count;
    (*Total) = min_count;
}
uint32_t calculate_hap_similarity(Hap_Align* p, uint32_t dir, uint32_t xCid, uint32_t yCid,
kvec_t_u32_warp* x_vecs, kvec_t_u32_warp* y_vecs, ma_hit_t_alloc* reverse_sources, 
asg_t *read_g, R_to_U* ruIndex, float Hap_rate, uint32_t seedOcc)
{
    uint32_t max_count = 0, min_count = 0;
    uint32_t xLen = x_vecs->a.n;
    uint32_t yLen = y_vecs->a.n;
    uint32_t xLeftBeg, xLeftLen, yLeftBeg, yLeftLen;
    uint32_t xRightBeg, xRightLen, yRightBeg, yRightLen;
    if(dir == 0)
    {
        xLeftBeg = 0; xLeftLen = p->q_pos; xRightBeg = p->q_pos; xRightLen = xLen - xRightBeg;
        yLeftBeg = 0; yLeftLen = p->t_pos; yRightBeg = p->t_pos; yRightLen = yLen - yRightBeg;
    }
    else
    {
        xLeftBeg = 0; xLeftLen = p->q_pos; xRightBeg = p->q_pos; xRightLen = xLen - xRightBeg;

        yLeftBeg = p->t_pos + 1; yLeftLen = yLen - yLeftBeg;
        yRightBeg = 0; yRightLen = p->t_pos + 1;
    }



    max_count = seedOcc;
    min_count = MIN(xLeftLen, yLeftLen) + MIN(xRightLen, yRightLen);
    if(min_count == 0) return NON_PLOID;
    if(max_count <= min_count*Hap_rate) return NON_PLOID;



    

    double xLeftMatch, xLeftTotal, yLeftMatch, yLeftTotal;
    double xRightMatch, xRightTotal, yRightMatch, yRightTotal;

    get_hap_similarity(x_vecs->a.a+xLeftBeg, xLeftLen, yCid, reverse_sources, read_g, ruIndex, 
    &xLeftMatch, &xLeftTotal);
    get_hap_similarity(y_vecs->a.a+yLeftBeg, yLeftLen, xCid, reverse_sources, read_g, ruIndex,
    &yLeftMatch, &yLeftTotal);

    get_hap_similarity(x_vecs->a.a+xRightBeg, xRightLen, yCid, reverse_sources, read_g, ruIndex, 
    &xRightMatch, &xRightTotal);
    get_hap_similarity(y_vecs->a.a+yRightBeg, yRightLen, xCid, reverse_sources, read_g, ruIndex,
    &yRightMatch, &yRightTotal);

    max_count = min_count = 0;
    if((xLeftMatch/xLeftTotal) >= (yLeftMatch/yLeftTotal))
    {
        max_count += xLeftMatch;
        min_count += xLeftTotal;
    }
    else
    {
        max_count += yLeftMatch;
        min_count += yLeftTotal;
    }

    if((xRightMatch/xRightTotal) >= (yRightMatch/yRightTotal))
    {
        max_count += xRightMatch;
        min_count += xRightTotal;
    }
    else
    {
        max_count += yRightMatch;
        min_count += yRightTotal;
    }

    if(min_count == 0) return NON_PLOID;
    if(max_count > min_count*Hap_rate) return PLOID;
	return NON_PLOID;
}


#define Hap_Align_Pos_key(a) ((((uint64_t)((a).q_pos))<<32)|((uint64_t)((a).t_pos)))
KRADIX_SORT_INIT(Hap_Align_Pos_sort, Hap_Align, Hap_Align_Pos_key, 8)

#define Hap_Align_Weight_key(a) ((a).is_color)
KRADIX_SORT_INIT(Hap_Align_Weight_sort, Hap_Align, Hap_Align_Weight_key, member_size(Hap_Align, is_color))

inline uint32_t merge_hap_hits(Hap_Align_warp* u_buffer)
{
    if(u_buffer->x.n == 0) return 0;

    radix_sort_Hap_Align_Pos_sort(u_buffer->x.a, u_buffer->x.a + u_buffer->x.n);
    uint32_t x_pos, x_pos_end, x_pos_beg;
    int k, i, j, m;

    /**
    for (i = 0; i < (int)u_buffer->x.n; i++)
    {
        fprintf(stderr, "****x: %u, y: %u, t_id: %u, is_color: %u\n", u_buffer->x.a[i].q_pos, u_buffer->x.a[i].t_pos,
        u_buffer->x.a[i].t_id, u_buffer->x.a[i].is_color);
    }
    **/
    
    i = u_buffer->x.n; i--;
    x_pos = u_buffer->x.a[i].q_pos; 
    x_pos_end = i;
    for (; i >= 0; i--)
    {
        k = i - 1;
        if(k < 0) continue;
        if(u_buffer->x.a[k].q_pos == u_buffer->x.a[i].q_pos &&
           u_buffer->x.a[k].t_pos+1 == u_buffer->x.a[i].t_pos)
        {
            u_buffer->x.a[k].is_color += u_buffer->x.a[i].is_color;
            u_buffer->x.a[i].t_id = (uint32_t)-1;
        }

        ///meet a new x_pos
        if(u_buffer->x.a[k].q_pos != u_buffer->x.a[i].q_pos)
        {
            x_pos_beg = i;
            for (m = x_pos_beg; m <= (int)x_pos_end; m++)
            {
                if(u_buffer->x.a[m].t_id == (uint32_t)-1) continue;
                for (j = k; j >= 0; j--)
                {
                    if(u_buffer->x.a[j].q_pos != x_pos - 1) break;
                    if(u_buffer->x.a[j].t_pos == u_buffer->x.a[m].t_pos ||
                       u_buffer->x.a[j].t_pos == u_buffer->x.a[m].t_pos - 1)
                    {
                        u_buffer->x.a[j].is_color += u_buffer->x.a[m].is_color;
                        u_buffer->x.a[m].t_id = (uint32_t)-1;
                        break;
                    }
                }
            }
            
            x_pos = u_buffer->x.a[k].q_pos;
            x_pos_end = k;
        }
    }


    

    for (i = 0, m = 0; i < (int)u_buffer->x.n; i++)
    {
        if(u_buffer->x.a[i].t_id!=(uint32_t)-1)
        {
            u_buffer->x.a[m].is_color = u_buffer->x.a[i].is_color;
            u_buffer->x.a[m].t_id = u_buffer->x.a[i].t_id;
            u_buffer->x.a[m].q_pos = u_buffer->x.a[i].q_pos;
            u_buffer->x.a[m].t_pos = u_buffer->x.a[i].t_pos;
            m++;
        }
    }

    u_buffer->x.n = m;

    /**
    for (i = 0; i < (int)u_buffer->x.n; i++)
    {
        fprintf(stderr, "####x: %u, y: %u, t_id: %u, is_color: %u\n", u_buffer->x.a[i].q_pos, u_buffer->x.a[i].t_pos,
        u_buffer->x.a[i].t_id, u_buffer->x.a[i].is_color);
    }
    **/
    radix_sort_Hap_Align_Weight_sort(u_buffer->x.a, u_buffer->x.a + u_buffer->x.n);
    return 0;
}

void get_hap_alignment(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources, 
buf_t* b_0, R_to_U* ruIndex, uint32_t* position_index, uint64_t* vote_counting, uint8_t* visit, 
kvec_t_u64_warp* u_vecs, Hap_Align_warp* u_buffer, kvec_t_u32_warp* x_vecs, kvec_t_u32_warp* y_vecs,
uint32_t cId, uint32_t contigNum, uint64_t* cBeg, float Hap_rate, asg_t *bi_g, uint32_t is_bi_edge)
{
    ma_utg_t* reads = NULL;
    uint32_t beg, end, i, j, k, rId, Hap_cId, y_cId, y_offset, qn, self_offset, uId, uOri, is_Unitig, xLen, seedOcc;
    uint64_t tmp;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    memset(vote_counting, 0, sizeof(uint64_t)*contigNum);
    memset(visit, 0, contigNum);
    u_vecs->a.n = x_vecs->a.n = y_vecs->a.n = 0;

    beg = (uint32_t)(cBeg[cId]);
    b_0->b.n = 0;
    get_unitig(ug->g, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);
    xLen = cBeg[cId]>>33;
    kv_resize(uint32_t, x_vecs->a, xLen);

    for (i = 0, self_offset = 0; i < b_0->b.n; i++)
    {
        uId = b_0->b.a[i]>>1;
        uOri = b_0->b.a[i]&(uint32_t)1;
        reads = &(ug->u.a[uId]);
        for (j = 0; j < reads->n; j++, self_offset++)
        {
            if(uOri == 1)
            {
                rId = reads->a[reads->n - j - 1]>>33;
            }
            else
            {
                rId = reads->a[j]>>33;
            }


            tmp = rId<<1;
            if(read_g->seq[rId].c == HAP_LABLE)
            {
                tmp = tmp | 1;
                kv_push(uint32_t, x_vecs->a, tmp);
                continue;
            }
            kv_push(uint32_t, x_vecs->a, tmp); 


            /**********************for debug*************************/
            ///debug_visit_flag(visit, contigNum);
            /**********************for debug*************************/

            qn = rId;
            for (k = 0; k < reverse_sources[qn].length; k++)
            {
                rId = Get_tn(reverse_sources[qn].buffer[k]);
            
                if(read_g->seq[rId].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                    if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
                }

                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_cId, &is_Unitig);
                if(is_Unitig == 0 || Hap_cId == (uint32_t)-1) continue;
                ///here rId is the id of the read coming from the different haplotype
                ///Hap_cId is the id of the corresponding contig (note here is the contig, instead of untig)
                if(visit[Hap_cId]!=0) continue;
                visit[Hap_cId] = 1;
                if(vote_counting[Hap_cId] < UINT64_MAX) vote_counting[Hap_cId]++;
            }

            reset_visit_flag(visit, read_g, ruIndex, contigNum, &(reverse_sources[qn]));
        }
    }


    u_vecs->a.n = 0;
    for (i = 0; i < contigNum; i++)
    {
        if(i == cId) continue;
        if(vote_counting[i] == 0) continue;
        tmp = vote_counting[i]; tmp = tmp << 32; tmp = tmp | (uint64_t)i;
        kv_push(uint64_t, u_vecs->a, tmp);
    }

    if(u_vecs->a.n == 0) return;

    sort_kvec_t_u64_warp(u_vecs, 1);


    ma_hit_t_alloc *x = NULL, *pre_x = NULL;
    Hap_Align* p = NULL;
    asg_arc_t *e = NULL;
    ///scan from the weightest candidate 
    for (i = 0; i < u_vecs->a.n; i++)
    {
        Hap_cId = (uint32_t)u_vecs->a.a[i];
        get_readSeq(ug, read_g, y_vecs, cBeg, b_0, Hap_cId);
        y_cId = Hap_cId;
        seedOcc = u_vecs->a.a[i]>>32;

        if(debug_purge_dup)
        {
            if(cId == 53)
            {
                fprintf(stderr, "cId: %u, y_cId: %u, seedOcc: %u\n", cId, y_cId, seedOcc);
            }
        }

        u_buffer->x.n = 0;
        for (k = 0; k < x_vecs->a.n; k++)
        {
            x = pre_x = NULL;
            rId = x_vecs->a.a[k]>>1;
            if(read_g->seq[rId].c == HAP_LABLE) continue;
            x = &(reverse_sources[rId]);

            if(k >= 1)
            {
                rId = x_vecs->a.a[k-1]>>1;
                if(read_g->seq[rId].c != HAP_LABLE)
                {
                    pre_x = &(reverse_sources[rId]);
                }
            }




            for (j = 0; j < x->length; j++)
            {
                rId = Get_tn(x->buffer[j]);
                if(read_g->seq[rId].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                    if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
                }

                ///there are two cases: 
                ///1. read at primary contigs, get_R_to_U() return its corresponding contig Id  
                ///2. read at alternative contigs, get_R_to_U() return (uint32_t)-1
                get_R_to_U(ruIndex, rId, &Hap_cId, &is_Unitig);
                if(is_Unitig == 0 || Hap_cId == (uint32_t)-1) continue;
                if(Hap_cId != y_cId) continue;

                y_offset = position_index[rId];
                ///if(retrieve_skip_overlaps(pre_x, rId) != (uint32_t)-1)
                if(retrieve_skip_overlaps(pre_x, Get_tn(x->buffer[j])) != (uint32_t)-1)
                {
                    check_duplicate(u_buffer, (uint32_t)-1, y_offset);
                    continue;
                } 
                if(check_duplicate(u_buffer, k, y_offset)) continue;

                kv_pushp(Hap_Align, u_buffer->x, &p);
                p->q_pos = k;
                p->t_pos = y_offset;
                p->t_id = y_cId;
                p->is_color = 1;
            }
        }

        merge_hap_hits(u_buffer);
        if(u_buffer->x.n == 0) continue;

        ///for (k = 0; k < u_buffer->x.n; k++)
        for (k = u_buffer->x.n-1; k >= 0; k--)
        {
            if(calculate_hap_similarity(&(u_buffer->x.a[k]), 0, cId, y_cId, x_vecs, y_vecs, 
            reverse_sources, read_g, ruIndex, Hap_rate, seedOcc)==PLOID)
            {
                break;
            }
            else if(calculate_hap_similarity(&(u_buffer->x.a[k]), 1, cId, y_cId, x_vecs, y_vecs, 
            reverse_sources, read_g, ruIndex, Hap_rate, seedOcc)==PLOID)
            {
                break;
            }
            if(k == 0)
            {
                k = (uint32_t)-1;
                break;
            }
        }

        ///if(k < u_buffer->x.n)
        if(k != (uint32_t)-1)
        {
            e = asg_arc_pushp(bi_g);
            e->del = 0;
            e->ol = 0;
            e->ul = cId; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(0);
            e->v = y_cId;

            if(is_bi_edge)
            {
                e = asg_arc_pushp(bi_g);
                e->del = 0;
                e->ol = 0;
                e->ul = y_cId; e->ul = e->ul << 32; e->ul = e->ul | (uint64_t)(0);
                e->v = cId;
            }
        }
    }
}

void print_gfa(asg_t *g)
{
    uint32_t v, i, n_vtx = g->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        if(g->seq[v>>1].del)
        {
            fprintf(stderr, "(D) v>>1: %u, v&1: %u, %.*s\n", v>>1, v&1, 
                (int)Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
            continue;
        } 

        fprintf(stderr, "(E) v>>1: %u, v&1: %u, %.*s\n", v>>1, v&1, 
                (int)Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1));
        
        uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            fprintf(stderr, "av[i].v: %u, av[i].ul: %u\n", 
                            av[i].v, (uint32_t)(av[i].ul>>32));
        }

    }
    
}



void print_purge_gfa(asg_t *g, uint64_t* cCount)
{
    uint32_t v, i, n_vtx = g->n_seq, beg, end;
    for (v = 0; v < n_vtx; v++)
    {
        if(g->seq[v>>1].del)
        {
            fprintf(stderr, "(D) v: %u, Len: %u, start>>1: %u, flag: %u\n", v, (uint32_t)(cCount[v]>>33),
            ((uint32_t)cCount[v])>>1, g->seq[v].len);
            continue;
        } 

        fprintf(stderr, "(E) v: %u, Len: %u, start>>1: %u, flag: %u\n", v, (uint32_t)(cCount[v]>>33),
        ((uint32_t)cCount[v])>>1, g->seq[v].len);
        
        uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            beg = av[i].ul>>32;
            end = av[i].v;

            fprintf(stderr, "****beg: %u (Len: %u) ---> end: %u (Len: %u)\n", 
            beg, (uint32_t)(cCount[beg]>>33),
            end, (uint32_t)(cCount[end]>>33));

        }

    }
    
}

void calculate_peak(uint64_t* counts, long long len, kvec_t_u32_warp* a, 
long long TotalMean, long long HapMean, long long* peak0, long long* peak1)
{
    if(len == 0) return;
    long long i, j, total_coverage = 0, current_coverage, step, pre, after, is_found = 0;
    for (i = 0; i < len; i++)
    {
        total_coverage = total_coverage + ((counts[i]>>32) * ((uint32_t)counts[i]));
    }
    current_coverage = total_coverage;
    for (i = len - 1; i >= 0; i--)
    {
        current_coverage = current_coverage - ((counts[i]>>32) * ((uint32_t)counts[i]));
        if(current_coverage <= total_coverage*0.9) break;
    }
    i++;
    step = i*0.1;
    if(step < 2) step = 2;
    ///fprintf(stderr, "step: %lld\n", step);
    

    a->a.n = 0;
    ///for (i = 0; i < len; i++)
    for (i = 1; i < len; i++)
    {
        pre = i - step/2;
        after = i + step/2;
        if(pre < 0) pre = 0;
        if(after >= len) after = len - 1;
        is_found = 0;
        for (j = pre; j < after; j++)
        {
            if(j == i) continue;
            if(((uint32_t)counts[i]) <= ((uint32_t)counts[j]))
            {
                is_found = 1;
                break;
            } 
        }
        
        if(is_found == 0)
        {
            kv_push(uint32_t, a->a, i);
        }
    }

    

    uint32_t peak_0_i, peak_1_i;
    peak_0_i = peak_1_i = 0;

    for (i = 0; i < (long long)a->a.n; i++)
    {
        if(((uint32_t)counts[a->a.a[peak_0_i]]) < ((uint32_t)counts[a->a.a[i]]))
        {
            peak_0_i = i;
        }
    }

    peak_1_i = (uint32_t)-1;
    for (i = 0; i < (long long)a->a.n; i++)
    {
        if(peak_0_i == i) continue;
        if(peak_1_i == (uint32_t)-1)
        {
            peak_1_i = i;
            continue;
        }
        if(((uint32_t)counts[a->a.a[peak_1_i]]) < ((uint32_t)counts[a->a.a[i]]))
        {
            peak_1_i = i;
        }
    }

    if(peak_0_i == peak_1_i) peak_1_i = (uint32_t)-1;

    
    // for (i = 0; i < (long long)a->a.n; i++)
    // {
    //     fprintf(stderr, "i:%lld, freq: %u, num: %u\n", i, 
    //     (uint32_t)(counts[a->a.a[i]]>>32), (uint32_t)counts[a->a.a[i]]);
    // }
    // fprintf(stderr, "peak_0_i: %u, peak_1_i: %u\n", peak_0_i, peak_1_i); 

    (*peak0) = (*peak1) = -1;
    /**
     * There are three cases:
     * 1. very low het rate: don't need purge_dup
     * 2. very high het rate: TotalMean should be equal to peak_0
     * 3. not such high het rate: two peaks, HapMean should be useful
     * **/
    
    if(peak_1_i == (uint32_t)-1)
    {
        (*peak0) = (uint32_t)(counts[a->a.a[peak_0_i]]>>32);
        return;
    }

    (*peak0) = (uint32_t)(counts[a->a.a[peak_0_i]]>>32);
    (*peak1) = (uint32_t)(counts[a->a.a[peak_1_i]]>>32);


    if(TotalMean >= (long long)((uint32_t)(counts[a->a.a[peak_0_i]]>>32)) && 
       TotalMean >= (long long)((uint32_t)(counts[a->a.a[peak_1_i]]>>32)))
    {
        (*peak1) = -1;
        return;
    }

    long long max_cov = MAX((*peak0), TotalMean);
    long long min_cov = MIN((*peak0), TotalMean);
    if(min_cov >= max_cov*0.8)
    {
        max_cov = MAX((*peak1), TotalMean);
        min_cov = MIN((*peak1), TotalMean);
        if(min_cov < max_cov*0.6)
        {
            (*peak1) = -1;
            return;
        }
    }


    long long hap_cov = ((uint32_t)(counts[a->a.a[MIN(peak_0_i, peak_1_i)]]>>32));
    max_cov = MAX(hap_cov, HapMean);
    min_cov = MIN(hap_cov, HapMean);
    ///fprintf(stderr, "hap_cov: %lld, max_cov: %lld, min_cov: %lld\n", hap_cov, max_cov, min_cov); 
    if(min_cov < max_cov*0.70)
    {
        (*peak1) = -1;
        return;
    }

    min_cov = MIN((*peak0), (*peak1));
    max_cov = MAX((*peak0), (*peak1));
    (*peak0) = min_cov;
    (*peak1) = max_cov;
}

void get_purge_coverage(asg_t *read_g, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut,
uint32_t* junk_cov, uint32_t* hap_cov, uint32_t* dip_cov)
{
    (*junk_cov) = (*hap_cov) = 0; (*dip_cov) = (uint32_t)-1;
    uint32_t i, j, num, m;
    long long R_bases = 0, C_bases = 0, T_R_bases = 0, T_C_bases = 0, total_coverage = 0, T_mean, Hap_mean;
    ma_hit_t *h = NULL;
    uint64_t* counts = (uint64_t*)calloc(read_g->n_seq, sizeof(uint64_t));
    kvec_t_u32_warp a;
    kv_init(a.a);

    for (i = 0; i < read_g->n_seq; ++i) 
    {
        R_bases = C_bases = 0;
        R_bases += coverage_cut[i].e - coverage_cut[i].s;
        for (j = 0; j < (uint64_t)(sources[i].length); j++)
        {
            h = &(sources[i].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
        if(R_bases == 0) continue;

        counts[i] = C_bases/R_bases;
        T_R_bases += R_bases;
        T_C_bases += C_bases;
        total_coverage += counts[i];
    }

    T_mean = 0;
    if(T_R_bases > 0) T_mean = T_C_bases/T_R_bases;

    radix_sort_arch64(counts, counts + read_g->n_seq);


    uint64_t tmp;
    i = m = 0;
    while(i < read_g->n_seq)
    {
        num = 0;
        for (j = i; j < read_g->n_seq; j++)
        {
            if(counts[i] != counts[j]) break;
            num++;
        }

        tmp = counts[i]; tmp = tmp << 32; tmp = tmp | (uint64_t)num;
        counts[m] = tmp;
        m++;
        i = j;
    }


    R_bases = C_bases = 0;
    for (i = 0; i < read_g->n_seq; ++i) 
    {
        if(read_g->seq[i].c != HAP_LABLE) continue;
        R_bases += coverage_cut[i].e - coverage_cut[i].s;
        for (j = 0; j < (uint64_t)(sources[i].length); j++)
        {
            h = &(sources[i].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }
    Hap_mean = 0;
    if(R_bases > 0) Hap_mean = C_bases/R_bases;



    long long peak0, peak1;
    calculate_peak(counts, m, &a, T_mean, Hap_mean, &peak0, &peak1);

   
    // fprintf(stderr, "T_mean: %lld, Hap_mean: %lld, peak0: %lld, peak1: %lld\n", 
    // T_mean, Hap_mean, peak0, peak1);


    (*junk_cov) = MIN(JUNK_COV, (peak0/2));
    if(peak1 == -1)
    {
        if(peak0 <= T_mean)
        {
            (*hap_cov) = peak0*1.5;
            (*dip_cov) = peak0*4.5;
        }
        else
        {
            (*hap_cov) = peak0*1.2;
            (*dip_cov) = peak0*3.6;
        }
    }
    else
    {
        num = MAX(T_mean, (peak0+peak1)/2);
        (*hap_cov) = num;
        (*dip_cov) = num*3;
    }
    
    free(counts);
    kv_destroy(a.a);
}


uint32_t get_single_coverage(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint32_t rId)
{
    long long R_bases = 0, C_bases = 0;
    uint32_t j;
    ma_hit_t *h = NULL;

    R_bases += coverage_cut[rId].e - coverage_cut[rId].s;
    for (j = 0; j < (uint64_t)(sources[rId].length); j++)
    {
        h = &(sources[rId].buffer[j]);
        C_bases += Get_qe((*h)) - Get_qs((*h));
    }

    if(R_bases == 0) return 0;
    return (C_bases/R_bases);
}


void further_clean_untig_graph_trio_advance(ma_ug_t *ug, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex, buf_t* b_0, 
uint8_t* visit, float density, uint32_t bi_graph_Len, uint32_t long_hap_overlap, float lable_match_rate)
{
    asg_t *bi_g = NULL;
    bi_g = asg_init();
    kvec_t(uint64_t) a;
    kv_init(a);
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, beg = 0, end, uId, uOri, cId, rId, i, j, k, self_offset, self_label_offset;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen;
    uint64_t uInfor = 0;
    ma_utg_t* reads;
    ///int flag;
    memset(visit, 0, nsg->n_seq);
    uint32_t* position_index = (uint32_t*)malloc(sizeof(uint32_t)*read_g->n_seq);
    memset(position_index, -1, sizeof(uint32_t)*read_g->n_seq);
    uint64_t* vote_counting = NULL;
    Hap_Align_warp u_buffer;
    kv_init(u_buffer.x);
    kvec_t_u32_warp x_vecs;
    kv_init(x_vecs.a);
    kvec_t_u32_warp y_vecs;
    kv_init(y_vecs.a);
    uint32_t junk_cov, hap_cov, dip_cov, junk_occ, repeat_occ, single_cov;


    /****************************may have bugs********************************/
    for (v = 0; v < nsg->n_seq; ++v) 
    {
        uId = v;
        if(nsg->seq[uId].c != HAP_LABLE) continue;
        reads = &(ug->u.a[uId]);
        for (k = 0; k < reads->n; k++)
        {
            rId = reads->a[k]>>33;
            read_g->seq[rId].c = HAP_LABLE;
        }
    }
    /****************************may have bugs********************************/

    
    get_purge_coverage(read_g, sources, coverage_cut, &junk_cov, &hap_cov, &dip_cov);
    fprintf(stderr, "junk_cov: %u, hap_cov: %u, dip_cov: %u\n", junk_cov, hap_cov, dip_cov);
    ///junk_cov = 50; dip_cov = 200;
    
    /******************************set ruIndex*********************************/
    for (v = 0; v < n_vtx; ++v) 
    {
        beg = v;
        if(nsg->seq[v>>1].c == ALTER_LABLE || nsg->seq[beg>>1].del || visit[beg>>1])
        {
            continue;
        }

        if(get_real_length(nsg, beg^1, NULL) == 1)
        {
            get_real_length(nsg, beg^1, &end);
            if(get_real_length(nsg, end^1, NULL) == 1) continue;
        }

        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);

        uInfor = 0;
        //scan all untigs
        for (i = 0; i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            visit[uId] = 1;
            reads = &(ug->u.a[uId]);
            uInfor += reads->n;
        }
        /****************************may have bugs********************************/
        ///uInfor = uInfor << 32; uInfor = uInfor | (uint64_t)beg;
        uInfor = uInfor << 33; uInfor = uInfor | (uint64_t)beg;
        /****************************may have bugs********************************/
        kv_push(uint64_t, a, uInfor);
    }

    ///sort by number of reads in a contig
    radix_sort_arch64(a.a, a.a + a.n);
    for (i = 0; i < (a.n>>1); ++i) 
    {
        uInfor = a.a[i];
        a.a[i] = a.a[a.n - i - 1];
        a.a[a.n - i - 1] = uInfor;
    }


    for (v = 0; v < a.n; v++)
    {
        ///all untig ID of this contig
        beg = (uint32_t)(a.a[v]);

        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);
       
        cId = v;
        ///individual contig
        if(get_real_length(nsg, beg^1, NULL) == 0 && get_real_length(nsg, end, NULL) == 0)
        {
            a.a[v] = a.a[v] | (uint64_t)(0x100000000);
        }
        ///set the contig Id for each read
        for (i = 0, self_offset = 0, self_label_offset = 0, junk_occ = 0, repeat_occ = 0; 
                                                                            i < b_0->b.n; i++)
        {
            uId = b_0->b.a[i]>>1;
            uOri = b_0->b.a[i]&(uint32_t)1;
            reads = &(ug->u.a[uId]);
            for (j = 0; j < reads->n; j++, self_offset++)
            {
                if(uOri == 1)
                {
                    rId = reads->a[reads->n - j - 1]>>33;
                }
                else
                {
                    rId = reads->a[j]>>33;
                }
                
                set_R_to_U(ruIndex, rId, cId, 1);
                position_index[rId] = self_offset;

                if(read_g->seq[rId].c == HAP_LABLE) self_label_offset++;
                single_cov = get_single_coverage(sources, coverage_cut, rId);
                if(single_cov <= junk_cov) junk_occ++;
                if(single_cov > dip_cov) repeat_occ++;

            }
        }


        if(junk_occ >= self_offset*DISCARD_RATE || repeat_occ >= self_offset*DISCARD_RATE)
        {
            for(i = 0; i < b_0->b.n; i++)
            {
                uId = b_0->b.a[i]>>1;
                uOri = b_0->b.a[i]&(uint32_t)1;
                reads = &(ug->u.a[uId]);
                for (j = 0; j < reads->n; j++, self_offset++)
                {
                    if(uOri == 1)
                    {
                        rId = reads->a[reads->n - j - 1]>>33;
                    }
                    else
                    {
                        rId = reads->a[j]>>33;
                    }

                    ruIndex->index[rId] = (uint32_t)-1;                    
                    position_index[rId] = (uint32_t)-1;
                }
            }

            asg_seq_set(bi_g, cId, BLACK, 0);
        }
        else if(self_offset >= long_hap_overlap && self_label_offset >= (self_offset*lable_match_rate))
        {
            asg_seq_set(bi_g, cId, RED, 0);
        }
        else
        {
            asg_seq_set(bi_g, cId, UNVISIT, 0);
        }
    }
    /******************************set ruIndex*********************************/

    



 
    vote_counting = (uint64_t*)malloc(sizeof(uint64_t)*a.n);
    memset(vote_counting, 0, sizeof(uint64_t)*a.n);
    ///scan all contigs from the longest one to the shortest one
    for (i = 0; i < a.n; i++)
    {
        if(debug_purge_dup)
        {
            if(i == 53)
            {
                fprintf(stderr, "+i: %u, Len: %u, start>>1: %u, flag: %u\n", i, 
                (uint32_t)(a.a[i]>>33), ((uint32_t)a.a[i])>>1, bi_g->seq[i].len);
            }
        }
        if(bi_g->seq[i].len == BLACK) continue;
        get_hap_alignment(ug, read_g, reverse_sources, b_0, ruIndex, position_index, 
        vote_counting, visit, &u_vecs, &u_buffer, &x_vecs, &y_vecs, i, a.n, a.a, 
        density, bi_g, 1);

        if(debug_purge_dup)
        {
            if(i == 53)
            {
                fprintf(stderr, "+i: %u, Len: %u, start>>1: %u, flag: %u\n", v, 
                (uint32_t)(a.a[i]>>33), ((uint32_t)a.a[i])>>1, bi_g->seq[i].len);
            }
        }

    }

    fprintf(stderr, "\n*****0*****\n");
    fprintf(stderr, "a.n: %u\n", (uint32_t)(a.n));
    process_bi_graph(bi_g);
    fprintf(stderr, "*****1*****\n");
    if(debug_purge_dup) print_purge_gfa(bi_g, a.a);
    fprintf(stderr, "*****2*****\n");
    bi_paration(bi_g, a.a, bi_graph_Len);
    fprintf(stderr, "*****3*****\n");

    for (v = 0; v < bi_g->n_seq; v++)
    {
        beg = (uint32_t)(a.a[v]);

        b_0->b.n = 0;
        get_unitig(nsg, ug, beg, &end, &nodeLen, &baseLen, &max_stop_nodeLen, &max_stop_baseLen, 1, b_0);

        if(bi_g->seq[v].len == BLACK)
        {
            ///if(debug_purge_dup) fprintf(stderr, "BLACK v: %u, beg>>1: %u\n", v, beg>>1);
            for (i = 0; i < b_0->b.n; i++)
            {
                nsg->seq[(b_0->b.a[i])>>1].c = ALTER_LABLE;
            }

            for (i = 0; i < b_0->b.n; i++)
            {
                asg_seq_drop(nsg, ((b_0->b.a[i])>>1));
            }
        }
        







        // fprintf(stderr, "cId: %u, beg>>1: %u, end>>1: %u, b_0->b.n: %u, Len: %u, type: %u\n", 
        // v, beg>>1, end>>1, (uint32_t)b_0->b.n, (uint32_t)(a.a[v]>>33), bi_g->seq[v].len);
        // uint32_t nv = asg_arc_n(bi_g, v), w;
        // asg_arc_t *av = asg_arc_a(bi_g, v);
        // for (i = 0; i < nv; ++i)
        // {
        //     w = av[i].v;
        //     fprintf(stderr, "w: %u\n", w);
        // }
    }

    
    /**
    for (i = 1; i < a.n; i++)
    {
        if((a.a[i]>>33) > (a.a[i-1]>>33))
        {
            fprintf(stderr, "hehe\n");
        }
    }

    
    for (i = 0; i < a.n; i++)
    {
        uint32_t k, get_cId, tLen = 0, is_Unitig;
        beg = (uint32_t)(a.a[i]);
        b_0->b.n = 0;ug->u.a[beg>>1].circ = 0;
        get_long_tip_length(nsg, &(ug->u), beg, &end, b_0);
        for (j = 0; j < b_0->b.n; j++)
        {
            uId = b_0->b.a[j]>>1;
            reads = &(ug->u.a[uId]);
            tLen += reads->n;
            for (k = 0; k < reads->n; k++)
            {
                rId = reads->a[k]>>33;
                get_R_to_U(ruIndex, rId, &get_cId, &is_Unitig);
                if(is_Unitig != 1 || get_cId != i)
                {
                    fprintf(stderr, "###is_Unitig: %u, get_cId: %u, i: %u\n", 
                    is_Unitig, get_cId, i);
                }
            }
        }

        uint32_t qLen = 0, m;
        for (m = 0; m < ruIndex->len; m++)
        {
            get_R_to_U(ruIndex, m, &get_cId, &is_Unitig);
            if(get_cId == (uint32_t)-1 || is_Unitig != 1) continue;
            if(get_cId == i)
            {
                qLen = 0;
                for (j = 0; j < b_0->b.n; j++)
                {
                    uId = b_0->b.a[j]>>1;
                    reads = &(ug->u.a[uId]);
                    tLen += reads->n;
                    for (k = 0; k < reads->n; k++)
                    {
                        rId = reads->a[k]>>33;
                        if(m == rId)
                        {
                            qLen = 1;
                            goto found;
                        }
                    }
                }
                found:;
                if(qLen == 0)
                {
                    fprintf(stderr, "***is_Unitig: %u, get_cId: %u, m: %u\n", 
                    is_Unitig, get_cId, m);
                }
            }
        }
    }
    **/

   uint32_t is_Unitig;
   for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    asg_cleanup(nsg);
    free(a.a);
    kv_destroy(u_vecs.a);
    kv_destroy(u_buffer.x);
    kv_destroy(x_vecs.a);
    kv_destroy(y_vecs.a);
    asg_destroy(bi_g);
    free(position_index);
    free(vote_counting);
}


void clean_untig_graph(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean)
{
    asg_t *g = ug->g;
    asg_cut_tip_primary(g, ug, tipsLen);
    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        ///need consider tangles
        asg_pop_bubble_primary(g, bubble_dist);
        ///need consider tangles
        asg_arc_cut_long_tip_primary(g, ug, tip_drop_ratio);
        ///need consider tangles
        ///note we need both the read graph and the untig graph
        untig_asg_arc_cut_long_equal_tips_assembly(ug, read_g, reverse_sources, 2, ruIndex);
        untig_asg_arc_cut_long_tip_primary_complex(ug, tip_drop_ratio, stops_threshold);
        untig_asg_arc_cut_long_equal_tips_assembly_complex(ug, read_g, reverse_sources, 2, 
        stops_threshold, ruIndex);
        if(is_final_clean)
        {
            untig_asg_arc_cut_chimeric(ug, read_g, reverse_sources, 2, stops_threshold, chimeric_rate,
            ruIndex);
        }
        cur_cons = get_graph_statistic(g);
    }

    asg_cut_tip_primary(g, ug, tipsLen);
    untig_asg_arc_simple_large_bubbles(ug, read_g, reverse_sources, 2, ruIndex);

    if(is_final_clean)
    {
        lable_hap_asg_by_ug(ug, read_g);
        further_clean_untig_graph(ug, read_g, reverse_sources, ruIndex, b_0, visit, density, miniHapLen,
        miniBiGraph, 200, 0.4, 0.5);
        adjust_asg_by_ug(ug, read_g);
    }
    
}


void clean_untig_graph_bubbles(ma_ug_t *ug, asg_t *read_g, ma_hit_t_alloc* reverse_sources,
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, buf_t* b_0, uint8_t* visit, float density, uint32_t miniHapLen, 
uint32_t miniBiGraph, float chimeric_rate, int is_final_clean)
{
    asg_t *g = ug->g;
    asg_cut_tip_primary(g, ug, tipsLen);
    long long pre_cons = get_graph_statistic(g);
    long long cur_cons = 0;
    while(pre_cons != cur_cons)
    {
        pre_cons = get_graph_statistic(g);
        ///need consider tangles
        asg_pop_bubble_primary(g, bubble_dist);
        cur_cons = get_graph_statistic(g);
    }

    asg_cut_tip_primary(g, ug, tipsLen);
    untig_asg_arc_simple_large_bubbles(ug, read_g, reverse_sources, 2, ruIndex);

    if(is_final_clean)
    {
        lable_hap_asg_by_ug(ug, read_g);
        adjust_asg_by_ug(ug, read_g);
    }
    
}



void resolve_simple_case(ma_ug_t *ug, asg_t* nsg)
{
    uint32_t v, n_vtx = nsg->n_seq * 2, nw, nv, w1, w2, i;
    asg_arc_t *av, *aw;
    ma_utg_t *v_x = NULL, *v_y = NULL;

    for (v = 0; v < n_vtx; ++v) 
    {
        if (nsg->seq[v>>1].del || nsg->seq[v>>1].c == ALTER_LABLE) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);

        ///for simple circle
        if(w1 == (w2^1))
        {
            ///fprintf(stderr, "* v>>1: %u\n", v>>1);
            EvaluateLen(ug->u, w1>>1) = EvaluateLen(ug->u, w1>>1) + EvaluateLen(ug->u, v>>1);
            ///nsg->seq[w1>>1].len = nsg->seq[w1>>1].len + nsg->seq[v>>1].len;

            v_x = &(ug->u.a[w1>>1]);
            v_y = &(ug->u.a[v>>1]);
            append_ma_utg_t(v_x, v_y);
            asg_seq_del(nsg, v>>1);
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 1 && get_real_length(nsg, w1, NULL) != 2) continue;

            ///fprintf(stderr, "# v>>1: %u\n", v>>1);
            EvaluateLen(ug->u, w1>>1) = EvaluateLen(ug->u, w1>>1) + EvaluateLen(ug->u, v>>1);
            ///nsg->seq[w1>>1].len = nsg->seq[w1>>1].len + nsg->seq[v>>1].len;
            
            v_x = &(ug->u.a[w1>>1]);
            v_y = &(ug->u.a[v>>1]);

            asg_seq_del(nsg, v>>1);
            if(get_real_length(nsg, w1, NULL) == 1) continue;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            av = asg_arc_a(nsg, w1^1);
            for (i = 0; i < nw; i++)
            {
                if(!aw[i].del)
                {
                    av[0].del = 0;aw[i].del = 1;
                    av[0].el = aw[i].el;
                    av[0].no_l_indel = aw[i].no_l_indel;
                    av[0].ol = aw[i].ol;
                    av[0].strong = aw[i].strong;
                    av[0].v = aw[i].v;
                    av[0].ul = (aw[i].ul)^(0x100000000);


                    av = asg_arc_a(nsg, aw[i].v^1);
                    nv = asg_arc_n(nsg, aw[i].v^1);
                    uint32_t k = 0;
                    for (k = 0; k < nv; k++)
                    {
                        if(av[k].v == (aw[i].ul>>32^1))
                        {
                            av[k].v = av[k].v^1;
                            break;
                        }
                    }
                    
                    break;
                }
            }

        }
    }
}

void print_node(asg_t* g, ma_ug_t *ug)
{
    int input_iv;
    uint32_t iv, v;
    asg_arc_t *av, *aw;
    uint32_t nv, nw, i, k, w;
    while (1)
    {
        fprintf(stderr, "\n\ninput v: ");
        if(scanf("%d", &input_iv) == 0) break;
        if(input_iv == -1) break;
        iv = input_iv;
        iv = iv << 1;

        v = iv;
        av = asg_arc_a(g, v);
        nv = asg_arc_n(g, v);
        fprintf(stderr, "0**************v>>1: %u, v&1: %u, c: %u, del: %u, len: %u**************\n", 
        v>>1, v&1, (uint32_t)g->seq[v>>1].c, (uint32_t)g->seq[v>>1].del, g->seq[v>>1].len);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = av[i].v;
            fprintf(stderr, "(%u) w>>1: %u, w&1: %u, ol: %u, eLen: %u\n", 
            i, w>>1, w&1, av[i].ol, (uint32_t)av[i].ul);
            

            aw = asg_arc_a(g, w^1);
            nw = asg_arc_n(g, w^1);
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == (v^1)) break;
            }

            fprintf(stderr, "self>>1: %u, self&1: %u, out>>1: %u, out^1: %u, is_sym: %u\n", 
            (uint32_t)(av[i].ul>>33), (uint32_t)((av[i].ul>>32)&1), av[i].v>>1, av[i].v&1, 
            (uint32_t)(k != nw));

        }






        v = iv^1;
        av = asg_arc_a(g, v);
        nv = asg_arc_n(g, v);
        fprintf(stderr, "1**************v>>1: %u, v&1: %u, c: %u, del: %u, len: %u**************\n", 
        v>>1, v&1, (uint32_t)g->seq[v>>1].c, (uint32_t)g->seq[v>>1].del, (uint32_t)g->seq[v>>1].len);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            w = av[i].v;
            fprintf(stderr, "w>>1: %u, w&1: %u, ol: %u, eLen: %u\n", 
            w>>1, w&1, (uint32_t)av[i].ol, (uint32_t)av[i].ul);
            

            aw = asg_arc_a(g, w^1);
            nw = asg_arc_n(g, w^1);
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == (v^1)) break;
            }

            fprintf(stderr, "self>>1: %u, self&1: %u, out>>1: %u, out^1: %u, is_sym: %u\n", 
            (uint32_t)(av[i].ul>>33), (uint32_t)((av[i].ul>>32)&1), av[i].v>>1, av[i].v&1, 
            (uint32_t)(k != nw));

        }

        if(ug != NULL)
        {
            ma_utg_t *ma_v = &(ug->u.a[v>>1]);
            fprintf(stderr, "\n###\nma_v->n: %u, ma_v->len: %u\n", ma_v->n, ma_v->len);
            fprintf(stderr, "start>>1: %u, start&1: %u, end>>1: %u, end&1: %u\n", 
            ma_v->start>>1, ma_v->start&1, ma_v->end>>1, ma_v->end&1);
            for (i = 0; i < ma_v->n; i++)
            {
                v = ma_v->a[i]>>32;
                fprintf(stderr, "i: %u, v>>1: %u, v&1: %u, len: %u\n", 
                i, v>>1, v&1, (uint32_t)ma_v->a[i]);
            }
            
        }
    }
}

void label_tangles(asg_t *sg, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, long long bubble_dist,
long long tipsLen, float tip_drop_ratio, long long stops_threshold, R_to_U* ruIndex, int just_bubble)
{
    double startTime = Get_T();
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));
    uint32_t v, sv, n_vtx, beg, end, next_uID = (uint32_t)-1;
    kvec_t_u32_warp u_vecs;
    kv_init(u_vecs.a);
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);

    ///for each untig, all node have the same direction
    ///and all node except the last one just have one edge
    ///the last one may have multiple edges
    ///for the useful untig, the signal is (u->n >= LongUntigThreshold && !u->circ)
    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
        IsMerge(ug->u, v) = 0;
    }

    resolve_simple_case(ug, nsg);
    asg_cleanup(nsg);
    asg_symm(nsg);

    ///print_node(nsg);

    if(just_bubble)
    {
        clean_untig_graph_bubbles(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, NULL, 0.8, 20, 200, 0.05, 0);
    }
    else
    {
        clean_untig_graph(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, NULL, 0.8, 20, 200, 0.05, 0);
    }
    
    

    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    uint32_t n_reduce, flag;
    n_vtx = nsg->n_seq * 2;
    while (1)
    {
        n_reduce = 0;
        for (v = 0; v < n_vtx; ++v) 
        {
            //as for return value: 0: do nothing, 1: unroll, 2: convex
            //we just need 1
            sv = v;
            flag = 0;
            while (1)
            {
                flag = walk_through(sg, ug, reverse_sources, minLongUntig, 
                    maxShortUntig, l_untig_rate, max_node_threshold, &b_0, &b_1, 
                    &u_vecs, visit, sv, &beg, &end, &next_uID, ruIndex);
                n_reduce += flag;
                if(flag != UNROLL_M)
                {
                    break;
                }
            }
        }
        if(n_reduce == 0) break;
    }

    asg_cleanup(nsg);
	asg_symm(nsg);

    if(just_bubble)
    {
        clean_untig_graph_bubbles(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, visit, 0.8, 20, 200, 0.05, 1);
    }
    else
    {
        clean_untig_graph(ug, sg, reverse_sources, bubble_dist, tipsLen, 
        tip_drop_ratio, stops_threshold, ruIndex, &b_0, visit, 0.8, 20, 200, 0.05, 1);
    }
    




    kv_destroy(u_vecs.a);
    ma_ug_destroy(ug);


    free(visit);
    free(b_0.b.a);
    free(b_1.b.a);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n", __func__, Get_T()-startTime);
    }
}

/*************************************for tangle resolve*************************************/

void recover_edges(asg_t* nsg, kvec_t_u64_warp* edges, uint64_t* nodes, uint64_t n,
uint32_t beg, uint32_t end, uint32_t beg_c, uint32_t end_c, uint32_t is_recover_edges)
{
    uint32_t i, v;
    asg_arc_t *a = NULL;
    for (i = 0; i < n; i++)
    {
        v = (uint32_t)nodes[i];
        nsg->seq[v].del = 0;
        nsg->seq[v].c = (nodes[i]>>32);
    }
    nsg->seq[beg>>1].c = beg_c; nsg->seq[beg>>1].del = 0;
    nsg->seq[end>>1].c = end_c; nsg->seq[end>>1].del = 0;

    if(is_recover_edges)
    {
        for (i = 0; i < edges->a.n; i++)
        {
            a = &(nsg->arc[(uint32_t)edges->a.a[i]]);
            a->del = 0;
        }
    }
}

///note: nodes does not have directions
void save_all_edges(kvec_t_u64_warp* u_vecs, asg_t* nsg, uint64_t* nodes, uint64_t n,
uint32_t beg, uint32_t end, uint32_t* beg_c, uint32_t* end_c)
{
    uint32_t i, v, nv, k;
    uint64_t tmp;
    asg_arc_t *av = NULL;
    u_vecs->a.n = 0;

    (*beg_c) = nsg->seq[beg>>1].c;
    (*end_c) = nsg->seq[end>>1].c;
    for (i = 0; i < n; i++)
    {
        v = (uint32_t)nodes[i];
        tmp = nsg->seq[v].c; tmp = tmp<<32; tmp = (tmp)|((uint64_t)v);
        nodes[i] = tmp;

        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;

        v = v<<1;
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            kv_push(uint64_t, u_vecs->a, (uint64_t)((uint64_t)av[k].ol << 32 | (av - nsg->arc + k)));   
        }
        

        v = v^1;
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            kv_push(uint64_t, u_vecs->a, (uint64_t)((uint64_t)av[k].ol << 32 | (av - nsg->arc + k)));   
        }
    }

    v = beg;
    if((!nsg->seq[v>>1].del)&&(nsg->seq[v>>1].c!=ALTER_LABLE))
    {
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==((uint32_t)nodes[i])) break;
            }
            if(i==n) continue;
            kv_push(uint64_t, u_vecs->a, (uint64_t)((uint64_t)av[k].ol << 32 | (av - nsg->arc + k)));   
        }
    }
    



    v = end^1;
    if((!nsg->seq[v>>1].del)&&(nsg->seq[v>>1].c!=ALTER_LABLE))
    {
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            for (i = 0; i < n; i++)
            {
                if((av[k].v>>1)==((uint32_t)nodes[i])) break;
            }
            if(i==n) continue;
            kv_push(uint64_t, u_vecs->a, (uint64_t)((uint64_t)av[k].ol << 32 | (av - nsg->arc + k)));   
        }
    }
}

///beg and end have directions, while nodes does not have
int drop_tips_at_tangle(asg_t* nsg, uint64_t* nodes, uint64_t n, uint32_t beg, uint32_t end, 
buf_t* bb)
{
    uint32_t i, v, k, m, convex, return_flag, next_tips = 1, reduce = 0;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    
    while(next_tips > 0)
    {
        next_tips = 0;

        for (i = 0; i < n; i++)
        {
            v = (uint32_t)nodes[i];
            if(nsg->seq[v].del) continue;
            if(nsg->seq[v].c == ALTER_LABLE) continue;
            v = v<<1;

            for (k = 0; k < 2; k++)
            {
                v = v|k;
                if(get_real_length(nsg, v^1, NULL) != 0) continue;

                bb->b.n = 0;
                return_flag = get_unitig(nsg, NULL, v, &convex, &ll, &tmp, &max_stop_nodeLen, 
                        &max_stop_baseLen, 1, bb);
                if(return_flag == LOOP) continue;
                if(return_flag == MUL_OUTPUT) next_tips++;

                for(m = 0; m < bb->b.n; m++)
                {
                    if(((bb->b.a[m]>>1)==(beg>>1)) || ((bb->b.a[m]>>1)==(end>>1))) break;
                }
                ///means this unitig consists of beg or end 
                if(m != bb->b.n) continue;

                for(m = 0; m < bb->b.n; m++)
                {
                    nsg->seq[(bb->b.a[m]>>1)].c = ALTER_LABLE;
                }

                for(m = 0; m < bb->b.n; m++)
                {
                    asg_seq_drop(nsg, (bb->b.a[m]>>1));
                }
                
                reduce++;
            }
        }
    }
    
    return reduce;
}

int pop_bubble_at_tangle(ma_ug_t *ug, int max_dist, 
uint64_t* nodes, uint64_t n, uint32_t beg, uint32_t end,
uint32_t positive_flag, uint32_t negative_flag)
{
    asg_t *g = ug->g;
    uint32_t i, v, k, n_vtx = g->n_seq*2;
    uint64_t n_pop = 0;
    buf_t b;
    if (!g->is_symm) asg_symm(g);
    memset(&b, 0, sizeof(buf_t));
    b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));

    v = beg;
    if((!g->seq[v>>1].del)&&(g->seq[v>>1].c!=ALTER_LABLE)&&get_real_length(g, v, NULL)>=2)
    {
        n_pop += asg_bub_pop1_primary_trio(ug->g, ug, v, max_dist, &b, positive_flag, negative_flag, 1);
    }

    v = end^1;
    if((!g->seq[v>>1].del)&&(g->seq[v>>1].c!=ALTER_LABLE)&&get_real_length(g, v, NULL)>=2)
    {
        n_pop += asg_bub_pop1_primary_trio(ug->g, ug, v, max_dist, &b, positive_flag, negative_flag, 1);
    }


    for (i = 0; i < n; i++)
    {
        v = (uint32_t)nodes[i];
        if(g->seq[v].del) continue;
        if(g->seq[v].c == ALTER_LABLE) continue;
        v = v<<1;

        for (k = 0; k < 2; k++)
        {
            v = v|k;
            if(get_real_length(g, v, NULL)<=1) continue;
            n_pop += asg_bub_pop1_primary_trio(ug->g, ug, v, max_dist, &b, positive_flag, negative_flag, 1);
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    ///if (n_pop) asg_cleanup(g);
    return n_pop;
}


int find_spec_node(asg_t *nsg, uint32_t vBeg, uint32_t vDest, uint32_t* forbidden_nodes, 
uint32_t forbidden_nodes_n, long long max_nodes, buf_t* bb, uint8_t* visit, ma_ug_t *ug)
{
    uint32_t nv, un_visit;
    asg_arc_t *av;

    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    uint32_t vEnd, i, k, v, is_found = 0, return_flag;
    long long nodeLen, baseLen, max_stop_nodeLen, max_stop_baseLen, totalNodes = 0;

    bb->b.n = 0;
    if(get_unitig(nsg, ug, vBeg, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb)==LOOP)
    {
        kdq_destroy(uint32_t, buf);

        for (i = 0; i < bb->b.n; i++)
        {
            if(bb->b.a[i] == vDest) break;
        }
        if(i != bb->b.n) return 1;

        return 0;
    }
    

    for (i = 0; i < bb->b.n; i++)
    {
        Set_vis(visit, bb->b.a[i], 0);
        ///Set_vis(visit, bb->b.a[i]^1, 0);
    }

    if(forbidden_nodes != NULL && forbidden_nodes_n != 0)
    {
        for (i = 0; i < forbidden_nodes_n; i++)
        {
            Set_vis(visit, forbidden_nodes[i], 0);
        }
    }
    
    

    v = vBeg;
    kdq_push(uint32_t, buf, v);

    while (kdq_size(buf) != 0)
    {

        v = *(kdq_pop(uint32_t, buf));

        bb->b.n = 0;
        return_flag = get_unitig(nsg, ug, v, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb);

        for (i = 0; i < bb->b.n; i++)
        {
            if(bb->b.a[i] == vDest) break;
        }

        if(i != bb->b.n) 
        {
            is_found = 1;
            break;
        }
        
        if(return_flag == LOOP) break;

        totalNodes += nodeLen;
        if(totalNodes > max_nodes) break;
        
        v = vEnd;
        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
        
            if(Get_vis(visit,av[k].v, 0)!=0) continue;

            un_visit = 0;
            bb->b.n = 0;
            get_unitig(nsg, ug, av[k].v, &vEnd, &nodeLen, &baseLen, 
                                        &max_stop_nodeLen, &max_stop_baseLen, 1, bb);
            for (i = 0; i < bb->b.n; i++)
            {
                if(Get_vis(visit, bb->b.a[i], 0) != 0) break;
            }
            if(i == bb->b.n) un_visit = 1;

            if(un_visit)
            {
                kdq_push(uint32_t, buf, av[k].v);
                for (i = 0; i < bb->b.n; i++)
                {
                    Set_vis(visit, bb->b.a[i], 0);
                }
            }
            
        }
    }

    kdq_destroy(uint32_t, buf);
    return is_found;
}

int drop_useless_edges(ma_ug_t *ug, buf_t* bb, uint8_t* visit, uint32_t beg, uint32_t end, 
uint64_t* nodes, uint64_t nodes_n, uint32_t NodesThres)
{
    asg_t *nsg = ug->g;
    asg_arc_t *av = NULL;
    uint32_t v, w, nv, i, k, dest, n_reduce = 0;
    uint32_t reach_beg[2];
    uint32_t reach_end[2];
    reach_beg[0] = reach_beg[1] = reach_end[0] = reach_end[1] = 0;

    v = beg; dest = end; 
    if(get_real_length(nsg, v, NULL) > 1)
    {
        w = v^1;
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(find_spec_node(nsg, av[k].v, dest, &w, 1, NodesThres, bb, visit, ug) == 0)
            {
                // if(get_real_length(nsg, v, NULL)<=1)
                // {
                //     fprintf(stderr, "false edge: beg: %u, end: %u\n", (uint32_t)(av[k].ul>>33), av[k].v>>1);
                // }
                
                av[k].del = 1;
                asg_arc_del(nsg, av[k].v^1, ((uint32_t)(av[k].ul>>32))^1, 1);
                n_reduce++;
            }
        }
    }
    


    v = end^1; dest = beg^1; 
    if(get_real_length(nsg, v, NULL) > 1)
    {
        w = v^1;
        av = asg_arc_a(nsg, v);
        nv = asg_arc_n(nsg, v);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(find_spec_node(nsg, av[k].v, dest, &w, 1, NodesThres, bb, visit, ug) == 0)
            {
                // if(get_real_length(nsg, v, NULL)<=1)
                // {
                //     fprintf(stderr, "false edge: beg: %u, end: %u\n", (uint32_t)(av[k].ul>>33), av[k].v>>1);
                // }

                av[k].del = 1;
                asg_arc_del(nsg, av[k].v^1, ((uint32_t)(av[k].ul>>32))^1, 1);
                n_reduce++;
            }
        }
    }

    uint32_t drop_dest, forbid_dest, drop_node;
    for (i = 0; i < nodes_n; i++)
    {
        v = (uint32_t)nodes[i];
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        v = v<<1;

        for (k = 0; k < 2; k++)
        {
            v = v|k;

            w = end;
            reach_beg[k] = find_spec_node(nsg, v, beg^1, &w, 1, NodesThres, bb, visit, ug);
            w = beg^1;
            reach_end[k] = find_spec_node(nsg, v, end, &w, 1, NodesThres, bb, visit, ug);
        }

        if((reach_beg[0]+reach_end[0])==0 || (reach_beg[1]+reach_end[1])==0) continue;

        k = drop_dest = forbid_dest = (uint32_t)-1;
        if((reach_beg[0]+reach_end[0])==1 && (reach_beg[1]+reach_end[1])>1) k = 0;
        if((reach_beg[0]+reach_end[0])>1 && (reach_beg[1]+reach_end[1])==1) k = 1;
        if(k == (uint32_t)-1) continue;

        drop_node = (uint32_t)nodes[i]; 
        drop_node = drop_node <<1; 
        drop_node = drop_node | (uint32_t)(1-k);
        if(reach_beg[k]==1)
        {
            //drop_dest = beg^1; forbid_dest = end;
            drop_dest = end; forbid_dest = beg^1;
        }

        if(reach_end[k]==1)
        {
            //drop_dest = end; forbid_dest = beg^1;
            drop_dest = beg^1; forbid_dest = end;
        }

        if(drop_dest == (uint32_t)-1 || forbid_dest == (uint32_t)-1) continue;

        av = asg_arc_a(nsg, drop_node);
        nv = asg_arc_n(nsg, drop_node);
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(find_spec_node(nsg, av[k].v, drop_dest, &forbid_dest, 1, NodesThres, bb, visit, ug) == 0)
            {
                av[k].del = 1;
                asg_arc_del(nsg, av[k].v^1, ((uint32_t)(av[k].ul>>32))^1, 1);

                // fprintf(stderr, "false inner edge: beg: %u, end: %u\n", 
                // (uint32_t)(av[k].ul>>33), av[k].v>>1);

                n_reduce++;
            }
        }

    }

    return n_reduce;
}

int detec_circles(asg_t *nsg, uint32_t vSrc, uint32_t vDest,  
uint8_t* visit, ma_ug_t *ug, uint64_t* nodes, uint64_t nodes_n)
{
    uint32_t nv, nw, knw, i, k, v, w;
    asg_arc_t *av, *aw;

    memset(visit, 0, nsg->n_seq);
    kdq_t(uint32_t) *buf;
    buf = kdq_init(uint32_t);

    vDest = vDest^1;

    Set_vis(visit, vSrc, 1);
    Set_vis(visit, vDest, 1);
    kdq_push(uint32_t, buf, vSrc);
    kdq_push(uint32_t, buf, vDest);

    for (i = 0; i < nodes_n; i++)
    {
        v = (uint32_t)nodes[i];
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        v = v<<1;

        for (k = 0; k < 2; k++)
        {
            v = v|k;
            if(get_real_length(nsg, v, NULL)==0)
            {
                Set_vis(visit, v, 1);
                kdq_push(uint32_t, buf, v^1);
            }
        }
    }



    while (kdq_size(buf) != 0)
    {
        v = *(kdq_pop(uint32_t, buf));
        ///if(Get_vis(visit, v, 1) == 0) fprintf(stderr, "ERROR\n");

        nv = asg_arc_n(nsg, v);
        av = asg_arc_a(nsg, v);
        for(k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            w = av[k].v^1;
            if(nsg->seq[w>>1].del) continue;
            if(Get_vis(visit, w, 1) != 0) continue;

            ///check indegree
            nw = asg_arc_n(nsg, w);
            aw = asg_arc_a(nsg, w);
            for (i = 0, knw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                if(Get_vis(visit, aw[i].v, 1) != 0) continue;
                knw++;
            }

            if(knw == 0)
            {
                kdq_push(uint32_t, buf, w^1);
                Set_vis(visit, (w^1), 1);
            }
        }
    }

    for (i = 0; i < nodes_n; i++)
    {
        v = (uint32_t)nodes[i];
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        v = v<<1;
        if(Get_vis(visit, v, 1) == 0) return 1;
    }

    return 0;
}

inline void check_connective(asg_t *nsg, uint32_t vBeg, uint32_t vEnd, long long totalNodeLen, buf_t* bb, 
uint8_t* visit, ma_ug_t *ug, uint64_t* nodes, uint64_t nodes_n, uint32_t* is_connect, 
uint32_t* is_circle)
{
    uint32_t w;
    (*is_connect) = (*is_circle) = 0;
    
    w = vBeg^1;
    (*is_connect) = find_spec_node(nsg, vBeg, vEnd, &w, 1, totalNodeLen, bb, visit, ug);
    (*is_circle) = detec_circles(nsg, vBeg, vEnd, visit, ug, nodes, nodes_n);
}


inline uint32_t process_tangles(ma_ug_t *ug, uint64_t* nodes, uint64_t nodes_n, uint32_t beg, uint32_t end,
buf_t* bb, uint8_t* visit, uint32_t trio_flag, long long totalNodeLen, long long totalBaseLen)
{
    uint32_t n_reduce = 1, v, w, kv, is_found = 0;
    asg_t* nsg = ug->g;

    while (n_reduce != 0)
    {
        while (n_reduce != 0)
        {
            n_reduce = 0;
            n_reduce += drop_tips_at_tangle(nsg, nodes, nodes_n, beg, end, bb);
            n_reduce += pop_bubble_at_tangle(ug, totalBaseLen, nodes, nodes_n, beg, end, trio_flag, DROP);
        }
        n_reduce = drop_useless_edges(ug, bb, visit, beg, end, nodes, nodes_n, totalNodeLen);
    }

    // if(pop_bubble_at_tangle(ug, 10000000, nodes, nodes_n, beg, end, trio_flag, DROP)!=0)
    // {
    //     fprintf(stderr, "false bubble popping: beg: %u, end: %u\n", beg>>1, end>>1);
    // }

    v = beg;
    is_found = 0;
    while (1)
    {
        if(v == end) is_found = 1;
        if(is_found == 1) break;

        kv = get_real_length(nsg, v, NULL);
        if(kv!=1) break;
        kv = get_real_length(nsg, v, &w);
        if(get_real_length(nsg, w^1, NULL)!=1) break;

        v = w;
        if(v == beg) break;
    }

    return is_found;
}

uint32_t cut_edges_progressive(ma_ug_t *ug, kvec_t_u64_warp* edges, float drop_ratio,
uint64_t* nodes, uint64_t nodes_n, uint32_t beg, uint32_t end, buf_t* bb, uint8_t* visit, 
uint32_t trio_flag, long long totalNodeLen, long long totalBaseLen)
{
    uint32_t k, v, i, kv, w, nv, ov_max, ban, is_found = 0;
    asg_arc_t *a = NULL, *av = NULL;
    asg_t* nsg = ug->g;
    radix_sort_arch64(edges->a.a, edges->a.a + edges->a.n);
    for (k = 0; k < edges->a.n; k++)
    {
        a = &nsg->arc[(uint32_t)edges->a.a[k]];
        if(a->del) continue;
        v = (a->ul)>>32; w = a->v;
        if(nsg->seq[v>>1].del || nsg->seq[v>>1].c == ALTER_LABLE) continue;
        if(nsg->seq[w>>1].del || nsg->seq[w>>1].c == ALTER_LABLE) continue;

        nv = asg_arc_n(nsg, v);
        if(nv<=1) continue;

        ov_max = 0;
        av = asg_arc_a(nsg, v);
        for(i = 0, kv = 0; i < nv; ++i) 
        {
            if(av[i].del) continue;
            if(ov_max < av[i].ol) ov_max = av[i].ol;
            ++kv;
        }

        if(kv<=1) continue;
        if(a->ol > ov_max * drop_ratio) continue;
        asg_arc_del(nsg, v, w, 1);
        asg_arc_del(nsg, w^1, v^1, 1);

        ban = beg^1;
        if(find_spec_node(nsg, beg, end, &ban, 1, totalNodeLen, bb, visit, ug) == 0)
        {
            is_found = 0;
            asg_arc_del(nsg, v, w, 0);
            asg_arc_del(nsg, w^1, v^1, 0);
            break;
        }

        is_found = process_tangles(ug, nodes, nodes_n, beg, end, bb, visit, trio_flag, totalNodeLen,
        totalBaseLen);
        if(is_found == 1) break;
    }

    return is_found;
}

///1. drop tips 2. break edges 3. drop tips 4. do bubble popping
///beg and end have direction, but nsu->a doesn't have
///note here we do everything on src, instead of ug
void unroll_tangle(ma_ug_t *ug, ma_ug_t *length_ug, uint64_t* nodes, uint64_t nodes_n, 
uint32_t beg, uint32_t end, kvec_t_u64_warp* edges, uint32_t trio_flag, buf_t* bb, uint8_t* visit,
float drop_ratio)
{
    uint32_t v, w, is_found = 0, convex, i, beg_c, end_c;
    long long tmp, max_stop_nodeLen, max_stop_baseLen, begLen, endLen, missLen;
    long long totalNodeLen = 0, totalBaseLen = 0;
    asg_t* nsg = ug->g;
    edges->a.n = 0;
    save_all_edges(edges, nsg, nodes, nodes_n, beg, end, &beg_c, &end_c);


    for (i = 0; i < nodes_n; i++)
    {
        v = (uint32_t)nodes[i];
        totalNodeLen += ug->u.a[v].n;
        totalBaseLen += ug->g->seq[v].len;
    }
    totalNodeLen += ug->u.a[beg>>1].n + ug->u.a[end>>1].n + 1;
    totalBaseLen += ug->g->seq[beg>>1].len + ug->g->seq[end>>1].len + 1;
    ///each might be visited twice
    totalNodeLen = totalNodeLen * 2;
    totalBaseLen = totalBaseLen * 2;

    /**************************debug**************************/
    // uint32_t debug_is_access = 0, debug_is_circle = 0;
    // check_connective(nsg, beg, end, totalNodeLen, bb, visit, ug, nodes, nodes_n, &debug_is_access, 
    // &debug_is_circle);
    // if(debug_is_access==0) fprintf(stderr, "Cannot approch end: beg: %u, end: %u\n", beg>>1, end>>1);
    // if(debug_is_circle==0) fprintf(stderr, "Not circle: beg: %u, end: %u\n", beg>>1, end>>1);
    // if(debug_is_circle==1) fprintf(stderr, "Circle: beg: %u, end: %u\n", beg>>1, end>>1);
    /**************************debug**************************/

    is_found = process_tangles(ug, nodes, nodes_n, beg, end, bb, visit, trio_flag, totalNodeLen, 
    totalBaseLen);

    ///if(is_found == 1) fprintf(stderr, "***Found: beg>>1: %u, end>>1: %u\n", beg>>1, end>>1);

    if(is_found == 0)
    {
        is_found = cut_edges_progressive(ug, edges, drop_ratio, nodes, nodes_n, 
        beg, end, bb, visit, trio_flag, totalNodeLen, totalBaseLen);

        ///if(is_found == 1) fprintf(stderr, "***Cutting Found: beg>>1: %u, end>>1: %u\n", beg>>1, end>>1);
    }
    

    if(is_found == 1)
    {
        get_unitig(length_ug->g, length_ug, beg^1, &convex, &begLen, &tmp, &max_stop_nodeLen, 
                        &max_stop_baseLen, 1, NULL);

        get_unitig(length_ug->g, length_ug, end, &convex, &endLen, &tmp, &max_stop_nodeLen, 
                        &max_stop_baseLen, 1, NULL);


        missLen = 0;
        for (i = 0; i < nodes_n; i++)
        {
            v = (uint32_t)nodes[i];
            if(nsg->seq[v].del || nsg->seq[v].c == ALTER_LABLE)
            {
                missLen += ug->u.a[v].n;
            }
        }
        ///resolved!
        if(missLen <= (begLen+endLen)*TANGLE_MISSED_THRES)
        {
            is_found = 0;
            w = beg^1;
            is_found = find_spec_node(nsg, beg, end, &w, 1, totalNodeLen, bb, visit, ug);
        }
        else
        {
            is_found = 0;
        }        
    }
    

    recover_edges(nsg, edges, nodes, nodes_n, beg, end, beg_c, end_c, 1-is_found);

}

void resolve_tangles(ma_ug_t *src, asg_t *read_g, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, R_to_U* ruIndex, uint32_t trio_flag, 
float drop_ratio)
{
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));
    uint32_t i, v, w, sv, n_vtx, beg, end, next_uID = (uint32_t)-1;
    kvec_t_u32_warp u_vecs;
    kv_init(u_vecs.a);

    kvec_t_u64_warp e_vecs;
    kv_init(e_vecs.a);


    
    ///note: we must reset start for each unitig
    n_vtx = src->g->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(src->g->seq[v].del) continue;
        if(src->u.a[v].m==0) continue;
        if(!(src->u.a[v].circ))
        {
            src->u.a[v].start = src->u.a[v].a[0]>>32;
        }
        else
        {
            src->u.a[v].start = UINT32_MAX;
        }
    }

    adjust_utg_advance(read_g, src, reverse_sources, ruIndex);
    ///note: we must reset start for each unitig
    n_vtx = src->g->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(src->g->seq[v].del) continue;
        if(src->u.a[v].m==0) continue;
        EvaluateLen(src->u, v) = src->u.a[v].n;
    }
    




    ma_ug_t *ug = NULL;
    ug = copy_untig_graph(src);

    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        ///nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
        IsMerge(ug->u, v) = 0;
    }

    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    uint32_t n_reduce, flag;
    n_vtx = nsg->n_seq * 2;
    
    while (1)
    {
        n_reduce = 0;
        for (v = 0; v < n_vtx; ++v) 
        {
            //as for return value: 0: do nothing, 1: unroll, 2: convex
            //we just need 1
            sv = v;
            flag = 0;
            while (1)
            {
                flag = walk_through(read_g, ug, reverse_sources, minLongUntig, 
                    maxShortUntig, l_untig_rate, max_node_threshold, &b_0, &b_1, 
                    &u_vecs, visit, sv, &beg, &end, &next_uID, ruIndex);
                n_reduce += flag;
                if(flag != UNROLL_M)
                {
                    break;
                }
            }
        }
        if(n_reduce == 0) break;
    }

    asg_cleanup(nsg);
	asg_symm(nsg);

    uint32_t uId, rId, is_Unitig, m, pre;
    ma_utg_t* nsu = NULL;
    for (v = 0; v < src->g->n_seq; v++)
    {
        uId = v;
        nsu = &(src->u.a[v]);
        if(nsu->m == 0) continue;
        if(src->g->seq[v].del) continue;
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            ///ori = nsu->a[k]>>32&1;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }


    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        nsu = &(ug->u.a[v]);
        m = 0;
        pre = (uint32_t)(-1);
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            get_R_to_U(ruIndex, rId, &uId, &is_Unitig);
            if(is_Unitig != 1 || uId == ((uint32_t)(-1))) continue;
            if(i > 0 && pre == uId) continue;
            pre = uId;
            nsu->a[m] = uId;
            m++;
        }
        nsu->n = m;
    }


    n_vtx = nsg->n_seq;
    for (i = 0; i < n_vtx; ++i) 
    {
        v = i;
        nsu = &(ug->u.a[v]);
        if(nsg->seq[v].del) continue;
        if(IsMerge(ug->u, v) == 0) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;

        v = v<<1;
        if(get_real_length(nsg, v, NULL) != 1) continue;
        get_real_length(nsg, v, &w);
        if(get_real_length(nsg, w^1, NULL) != 1) continue;
        if(IsMerge(ug->u, w>>1) != 0) continue;
        beg = w^1;


        v = v^1;
        if(get_real_length(nsg, v, NULL) != 1) continue;
        get_real_length(nsg, v, &w);
        if(get_real_length(nsg, w^1, NULL) != 1) continue;
        if(IsMerge(ug->u, w>>1) != 0) continue;
        end = w;


        
        

        ///we have three types of merged nodes
        ///1) CONVEX_M: one direction has two out-nodes, another direction has one out-node
        ///2) UNROLL_E: one direction has one out-node, another direction doesn't has out-node
        ///3) UNROLL_M: both directions have one out-node
        ///here we just need UNROLL_M
        ///beg and end must be unchanged in src
        unroll_tangle(src, ug, nsu->a, nsu->n, beg, end, &e_vecs, trio_flag, &b_0, visit, drop_ratio);

        /**
        v = v<<1;
        if(get_real_length(nsg, v, NULL) != 1)
        {
            fprintf(stderr, "\nERROR 1, v>>1: %u, get_real_length(nsg, v, NULL): %u\n", 
            v>>1, get_real_length(nsg, v, NULL));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        get_real_length(nsg, v, &w);
        if(get_real_length(nsg, w^1, NULL) != 1)
        {
            fprintf(stderr, "\nERROR 2, v>>1: %u, v&1: %u, get_real_length(nsg, w^1, NULL): %u\n", 
            v>>1, v&1, get_real_length(nsg, w^1, NULL));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        if(IsMerge(ug->u, w>>1) != 0)
        {
            fprintf(stderr, "\nERROR 3, v>>1: %u, v&1: %u, IsMerge(ug->u, w>>1): %u\n", v>>1, v&1, 
            IsMerge(ug->u, w>>1));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        beg = w^1;


        v = v^1;
        if(get_real_length(nsg, v, NULL) != 1)
        {
            fprintf(stderr, "\nERROR 4, v>>1: %u, v&1: %u, get_real_length(nsg, v, NULL): %u\n", 
            v>>1, v&1, get_real_length(nsg, v, NULL));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        get_real_length(nsg, v, &w);
        if(get_real_length(nsg, w^1, NULL) != 1)
        {
            fprintf(stderr, "\nERROR 5, v>>1: %u, v&1: %u, get_real_length(nsg, w^1, NULL): %u\n", 
            v>>1, v&1, get_real_length(nsg, w^1, NULL));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        if(IsMerge(ug->u, w>>1) != 0)
        {
            fprintf(stderr, "\nERROR 6, v>>1: %u, v&1: %u, IsMerge(ug->u, w>>1): %u\n", v>>1, v&1, 
            IsMerge(ug->u, w>>1));
            for (m = 0; m < nsu->n; m++)
            {
                fprintf(stderr, "uId: %lu, ", nsu->a[m]);
            }
            fprintf(stderr, "\n");
            continue;
        } 
        end = w;

        
        fprintf(stderr, "\n***\nbeg>>1: %u, beg&1: %u\n", beg>>1, beg&1);
        fprintf(stderr, "end>>1: %u, end&1: %u\n", end>>1, end&1);
        **/
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    kv_destroy(u_vecs.a);
    kv_destroy(e_vecs.a);
    
    ma_ug_destroy(ug);
    free(visit);
    free(b_0.b.a);
    free(b_1.b.a);

    
    ///note: we must reset start for each unitig
    n_vtx = src->g->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(src->g->seq[v].del) continue;
        if(src->u.a[v].m==0) continue;
        if(!(src->u.a[v].circ))
        {
            src->u.a[v].start = src->u.a[v].a[0]>>32;
        }
        else
        {
            src->u.a[v].start = UINT32_MAX;
        }
    }
    adjust_utg_advance(read_g, src, reverse_sources, ruIndex);
    ///note: we must reset start for each unitig
    n_vtx = src->g->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(src->g->seq[v].del) continue;
        if(src->u.a[v].m==0) continue;
        EvaluateLen(src->u, v) = src->u.a[v].n;
    }
    
}





void deduplicate(ma_ug_t *src, asg_t *read_g, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, R_to_U* ruIndex, uint32_t resolve_tangle)
{
    uint32_t i, v, sv, n_vtx, beg, end, next_uID = (uint32_t)-1, uId, is_Unitig, rId;
    ma_utg_t* nsu = NULL;
    ma_ug_t *ug = NULL;
    kvec_t_u32_warp u_vecs;
    kv_init(u_vecs.a);
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    ug = copy_untig_graph(src);

    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
        IsMerge(ug->u, v) = 0;
    }

    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    uint32_t n_reduce, flag;
    n_vtx = nsg->n_seq * 2;

    if(resolve_tangle)
    {
        while (1)
        {
            n_reduce = 0;
            for (v = 0; v < n_vtx; ++v) 
            {
                //as for return value: 0: do nothing, 1: unroll, 2: convex
                //we just need 1
                sv = v;
                flag = 0;
                while (1)
                {
                    flag = walk_through(read_g, ug, reverse_sources, minLongUntig, 
                        maxShortUntig, l_untig_rate, max_node_threshold, &b_0, &b_1, 
                        &u_vecs, visit, sv, &beg, &end, &next_uID, ruIndex);
                    n_reduce += flag;
                    if(flag != UNROLL_M)
                    {
                        break;
                    }
                }
            }
            if(n_reduce == 0) break;
        }

        asg_cleanup(nsg);
        asg_symm(nsg);
    }
    

    further_clean_untig_graph_trio(ug, read_g, reverse_sources, ruIndex, &b_0, visit, 
    0.8, 200, 200, 0.4, 0.5);

    for (v = 0; v < src->g->n_seq; v++)
    {
        uId = v;
        nsu = &(src->u.a[v]);
        if(nsu->m == 0) continue;
        if(src->g->seq[v].del) continue;
        if(src->g->seq[v].c==ALTER_LABLE) continue;
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }



    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c!=ALTER_LABLE) continue;
        nsu = &(ug->u.a[v]);
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            get_R_to_U(ruIndex, rId, &uId, &is_Unitig);
            if(is_Unitig != 1 || uId == ((uint32_t)(-1))) continue;
            src->g->seq[uId].c = ALTER_LABLE;
        }
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    kv_destroy(u_vecs.a);    
    ma_ug_destroy(ug);
    free(visit);
    free(b_0.b.a);
    free(b_1.b.a);
}


void deduplicate_advance(ma_ug_t *src, asg_t *read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long minLongUntig, 
long long maxShortUntig, float l_untig_rate, float max_node_threshold, R_to_U* ruIndex, 
uint32_t resolve_tangle)
{
    uint32_t i, v, sv, n_vtx, beg, end, next_uID = (uint32_t)-1, uId, is_Unitig, rId;
    ma_utg_t* nsu = NULL;
    ma_ug_t *ug = NULL;
    kvec_t_u32_warp u_vecs;
    kv_init(u_vecs.a);
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));

    ug = copy_untig_graph(src);

    asg_t* nsg = ug->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
        IsMerge(ug->u, v) = 0;
    }

    uint8_t* visit = NULL;
    visit = (uint8_t*)malloc(sizeof(uint8_t) * nsg->n_seq);
    uint32_t n_reduce, flag;
    n_vtx = nsg->n_seq * 2;

    if(resolve_tangle)
    {
        while (1)
        {
            n_reduce = 0;
            for (v = 0; v < n_vtx; ++v) 
            {
                //as for return value: 0: do nothing, 1: unroll, 2: convex
                //we just need 1
                sv = v;
                flag = 0;
                while (1)
                {
                    flag = walk_through(read_g, ug, reverse_sources, minLongUntig, 
                        maxShortUntig, l_untig_rate, max_node_threshold, &b_0, &b_1, 
                        &u_vecs, visit, sv, &beg, &end, &next_uID, ruIndex);
                    n_reduce += flag;
                    if(flag != UNROLL_M)
                    {
                        break;
                    }
                }
            }
            if(n_reduce == 0) break;
        }

        asg_cleanup(nsg);
        asg_symm(nsg);
    }
    

    further_clean_untig_graph_trio_advance(ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, &b_0, visit, 
    0.7, /**200, 200,**/50, 50, 0.5);

    for (v = 0; v < src->g->n_seq; v++)
    {
        uId = v;
        nsu = &(src->u.a[v]);
        if(nsu->m == 0) continue;
        if(src->g->seq[v].del) continue;
        if(src->g->seq[v].c==ALTER_LABLE) continue;
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }



    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c!=ALTER_LABLE) continue;
        nsu = &(ug->u.a[v]);
        for (i = 0; i < nsu->n; i++)
        {
            rId = nsu->a[i]>>33;
            get_R_to_U(ruIndex, rId, &uId, &is_Unitig);
            if(is_Unitig != 1 || uId == ((uint32_t)(-1))) continue;
            src->g->seq[uId].c = ALTER_LABLE;
        }
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    kv_destroy(u_vecs.a);    
    ma_ug_destroy(ug);
    free(visit);
    free(b_0.b.a);
    free(b_1.b.a);
}



/*************************************for tangle resolve*************************************/

uint32_t copy_ug_node(ma_ug_t *ug, asg_t* nsg, uint32_t v)
{
    ma_utg_t *p;
    ma_utg_t *o = &(ug->u.a[v]);
    uint32_t n_node = nsg->n_seq;
    asg_seq_set(nsg, n_node, nsg->seq[v].len, 0);
    kv_pushp(ma_utg_t, ug->u, &p);

    p->s = 0, p->start = o->start, p->end = o->end, p->len = o->len;
    p->n = o->n, p->circ = o->circ, p->m = o->m;
    p->a = (uint64_t*)malloc(8 * p->m);
    memcpy(p->a, o->a, (8*p->m));

    return n_node;
}


///v and w here have directions
uint32_t collect_ma_utg_ts(ma_ug_t *ug, uint32_t v, uint32_t w, ma_utg_t* result)
{
    uint32_t des_dir = w&1;
    long long i;
    uint64_t* aim = NULL;
    ma_utg_t *ma_w = &(ug->u.a[w>>1]);

    if(v == (uint32_t)-1)
    {
        result->start = ma_w->start;
        result->circ = ma_w->circ;
        result->end = ma_w->end;
        result->len = ma_w->len;
        result->n = 0;
    }

    
    if(result->m < (result->n + ma_w->n))
    {
        result->m = (result->n + ma_w->n);
        result->a = (uint64_t*)realloc(result->a, result->m * sizeof(uint64_t));
    }

    aim = result->a + result->n;

    if(des_dir == 1)
    {
        for (i = 0; i < (long long)ma_w->n; i++)
        {
            aim[ma_w->n - i - 1] = (ma_w->a[i])^(uint64_t)(0x100000000);
        }
    }
    else
    {
        for (i = 0; i < (long long)ma_w->n; i++)
        {
            aim[i] = ma_w->a[i];
        }
    }

    result->n = result->n + ma_w->n;
    return 1;
}


void debug_utg_graph(ma_ug_t *ug, asg_t* read_g, int require_equal_nv, int test_tangle)
{
    asg_t* nsg = ug->g;
    uint32_t n_vtx = nsg->n_seq, i, j, k, l, totalLen, v, nv, nw, w, untig_v, rid_v;
    asg_arc_t *aw = NULL, *av = NULL;
    for (i = 0; i < n_vtx; i++)
    {
        if(ug->g->seq[i].del) continue;

        totalLen = 0;
        ma_utg_t* result = &(ug->u.a[i]);
        if(result->n == 0) continue;
        v = (uint64_t)(result->a[0])>>32;
        if(result->start != UINT32_MAX && result->start != v) fprintf(stderr, "hehe\n");
        v = (uint64_t)(result->a[result->n-1])>>32;
        if(result->end != UINT32_MAX && result->end != (v^1)) fprintf(stderr, "haha\n");
        uint64_t end_index = result->n-1;
        if(result->start == UINT32_MAX && result->end == UINT32_MAX) end_index++;
        for (j = 0; j < end_index; j++)
        {
            v = (uint64_t)(result->a[j])>>32;
            if(j == result->n-1)
            {
                w = (uint64_t)(result->a[0])>>32; 
            }
            else
            {
                w = (uint64_t)(result->a[j+1])>>32; 
            }
            
            

            av = asg_arc_a(read_g, v);
            nv = asg_arc_n(read_g, v);
            l = 0;
            for (k = 0; k < nv; k++)
            {
                if(av[k].del) continue;
                if(av[k].v == w) 
                {
                    l = asg_arc_len(av[k]);
                    break;
                }
            }
            if(k == nv) fprintf(stderr ,"******error, j: %u, k: %u, nv: %u\n", j, k, nv);
            if(l != (uint32_t)(result->a[j])) fprintf(stderr ,"[E::%s]ERROR Length\n", __func__);
            totalLen = totalLen + l;
        }


        if(j < result->n)
        {
            v = (uint64_t)(result->a[j])>>32;
            l = read_g->seq[v>>1].len;
            if(l != (uint32_t)(result->a[j])) fprintf(stderr ,"*** ERROR Length, i: %u\n", i);
            totalLen = totalLen + l;
        }

        if(totalLen != result->len)
        {
            fprintf(stderr ,"[E::%s] ERROR Total Length, i: %u\n", __func__, i);
        } 
                


        if(ug->u.a[i].start == UINT32_MAX && ug->u.a[i].end == UINT32_MAX) continue;

        v = i<<1; v=v^1;
        av = asg_arc_a(ug->g, v);
        nv = asg_arc_n(ug->g, v);

        w = (ug->u.a[v>>1].start^1);
        aw = asg_arc_a(read_g, w);
        nw = asg_arc_n(read_g, w);


        if(require_equal_nv && get_real_length(ug->g, v, NULL) != get_real_length(read_g, w, NULL))
        {
            fprintf(stderr, "#########ERROR: i: %u, nv: %u, nw: %u\n", i, nv, nw);
        }

        for (j = 0; j < nv; j++)
        {
            if(av[j].del) continue;
            untig_v = av[j].v;
            if(untig_v&1) rid_v = ug->u.a[untig_v>>1].end;
            else rid_v = ug->u.a[untig_v>>1].start;
            
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == rid_v) break;
            }

            if(k == nw) fprintf(stderr, "#########ERROR: i: %u\n", i);
            if((k != nw) && (av[j].ol != aw[k].ol))
            {
                fprintf(stderr, "#########????????ERROR\n");
                fprintf(stderr, "nv: %u, nw: %u\n", nv, nw);
                fprintf(stderr, "av[%u].ol: %u, aw[%u].ol: %u, untig_v>>1: %u, untig_v&1: %u\n", 
                j, av[j].ol, k, aw[k].ol, untig_v>>1, untig_v&1);
            }
        }






        v = v^1;
        av = asg_arc_a(ug->g, v);
        nv = asg_arc_n(ug->g, v);

        w = (ug->u.a[v>>1].end^1);
        aw = asg_arc_a(read_g, w);
        nw = asg_arc_n(read_g, w);
        if(require_equal_nv && get_real_length(ug->g, v, NULL) != get_real_length(read_g, w, NULL))
        {
            fprintf(stderr, "*******ERROR: i: %u, nv: %u, nw: %u\n", i, nv, nw);
        }
        for (j = 0; j < nv; j++)
        {
            if(av[j].del) continue;
            untig_v = av[j].v;
            if(untig_v&1) rid_v = ug->u.a[untig_v>>1].end;
            else rid_v = ug->u.a[untig_v>>1].start;
            
            for (k = 0; k < nw; k++)
            {
                if(aw[k].del) continue;
                if(aw[k].v == rid_v) break;
            }

            if(k == nw) fprintf(stderr, "***********ERROR: i: %u\n", i);
            if((k != nw) && (av[j].ol != aw[k].ol))
            {
                fprintf(stderr, "***********????????ERROR\n");
                fprintf(stderr, "nv: %u, nw: %u\n", nv, nw);
                fprintf(stderr, "av[%u].ol: %u, aw[%u].ol: %u, untig_v>>1: %u, untig_v&1: %u\n", 
                j, av[j].ol, k, aw[k].ol, untig_v>>1, untig_v&1);
            }
        }

    }


    if(test_tangle != 1) return;
    fprintf(stderr, "test_tangle: %u\n", test_tangle);
    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t w1, w2, beg, end, rnw;
        if(nsg->seq[v>>1].del) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);

        ///for simple circle
        if(w1 == (w2^1))
        {
            beg = end = 0;
            aw = asg_arc_a(nsg, w1^1);
            nw = asg_arc_n(nsg, w1^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == (v^1)) continue;
                beg = aw[i].v;
            }
            if(rnw != 2) continue;

            aw = asg_arc_a(nsg, w2^1);
            nw = asg_arc_n(nsg, w2^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == v) continue;
                end = aw[i].v;
            }
            if(rnw != 2) continue;

            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            if((beg>>1) == (end>>1)) continue;
            fprintf(stderr, "\n************\n");
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 2) continue;
            end = beg = (uint32_t)-1;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            for (i = 0, rnw = 0; i < nw && rnw < 2; i++)
            {
                if(aw[i].del) continue;
                if(rnw == 0) beg = aw[i].v;
                if(rnw == 1) end = aw[i].v;
                rnw++;
            }
            if((beg>>1) == (end>>1)) continue;
            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            fprintf(stderr, "\n#################\n");
        }
    }
}

///just merge, don't delete anything
void merge_ug_nodes(ma_ug_t *ug, asg_t* read_g, kvec_t_u64_warp* array)
{
    if(array->a.n == 0) return;
    uint32_t i, k, v, w;
    ma_utg_t result;
    memset(&result, 0, sizeof(ma_utg_t));
    v = w = (uint32_t)-1;
    for (i = 0; i < array->a.n; i++)
    {
        w = array->a.a[i];
        collect_ma_utg_ts(ug, v, w, &result);
        v = w;
    }

    

    if(result.n == 0) return;
    asg_arc_t *av = NULL;
    uint32_t nv, l;
    result.len = 0;
    for (i = 0; i < result.n - 1; i++)
    {
        
        v = (uint64_t)(result.a[i])>>32;
        w = (uint64_t)(result.a[i + 1])>>32; 
        av = asg_arc_a(read_g, v);
        nv = asg_arc_n(read_g, v);
        l = 0;
        
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == w) 
            {
                l = asg_arc_len(av[k]);
                break;
            }
        }
        if(k == nv) fprintf(stderr ,"******error, i: %u, k: %u, nv: %u\n", i, k, nv);
        result.a[i] = v; result.a[i] = result.a[i]<<32; result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;
    }
    

    if(i < result.n)
    {
        v = (uint64_t)(result.a[i])>>32;
        l = read_g->seq[v>>1].len;
        result.a[i] = v;
        result.a[i] = result.a[i]<<32;
        result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;
    }
    //has already set result.a, result.len, result.n, result.m
    result.circ = 0;
    result.start = result.a[0]>>32;
    result.end = (result.a[result.n-1]>>32)^1;

    
    uint32_t beg_uid = array->a.a[0];
    uint32_t end_uid = array->a.a[array->a.n-1];
    uint32_t realLen = array->a.n;
    uint32_t new_uid = array->a.a[0]>>1;
    uint64_t kmp;





    /*******************************just for debug**********************************/
    /**
    if(beg_uid&1)
    {
        if(result.start != ug->u.a[beg_uid>>1].end)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    else
    {
        if(result.start != ug->u.a[beg_uid>>1].start)
        {
            fprintf(stderr, "ERROR\n");
        }
    }

    if(end_uid&1)
    {
        if(result.end != ug->u.a[end_uid>>1].start)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    else
    {
        if(result.end != ug->u.a[end_uid>>1].end)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    **/
    /*******************************just for debug**********************************/

    asg_arc_t *aw = NULL;
    uint32_t nw = 0;
    ///corresponding to direction 1 of new node
    v = beg_uid^1;
    av = asg_arc_a(ug->g, v);
    nv = asg_arc_n(ug->g, v);
    ///fprintf(stderr, "beg_uid_v>>1: %u, beg_uid_v&1: %u, nv: %u\n", v>>1, v&1, nv);
    for (k = 0; k < nv; k++)
    {
        if(av[k].del) continue;
        kmp = new_uid<<1; kmp = kmp^1; kmp = kmp << 32;

        ///if((av[k].v>>1) == (end_uid>>1)) continue;
        w = av[k].v^1; aw = asg_arc_a(ug->g, w); nw = asg_arc_n(ug->g, w);
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) break;
        }
        if(i == nw) fprintf(stderr, "[E::%s] ERROR at %s:%d\n", __func__, __FILE__, __LINE__);
        kmp = kmp | (uint64_t)(aw[i].ol);


        ///kmp = kmp | (uint64_t)(((ug->g)->idx[v]>>32) + k);/**kmp = kmp | av[k].v;**/
        ///here kmp is ul
        kv_push(uint64_t, array->a, kmp);
        kmp = av[k].ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(av[k].v);
        ///here kmp is ol + v
        kv_push(uint64_t, array->a, kmp);
        ///fprintf(stderr, "*av[%u].v>>1: %u, v&1: %u, ol: %u\n", k, av[k].v>>1, av[k].v&1, av[k].ol);
    }

    ///corresponding to direction 0 of new node
    v = end_uid;
    av = asg_arc_a(ug->g, v);
    nv = asg_arc_n(ug->g, v);
    ///fprintf(stderr, "end_uid_v>>1: %u, end_uid_v&1: %u, nv: %u\n", v>>1, v&1, nv);
    for (k = 0; k < nv; k++)
    {
        if(av[k].del) continue;
        kmp = new_uid<<1; kmp = kmp << 32;

        w = av[k].v^1; 
        aw = asg_arc_a(ug->g, w); 
        nw = asg_arc_n(ug->g, w);
        for (i = 0; i < nw; i++)
        {
            if(aw[i].del) continue;
            if(aw[i].v == (v^1)) break;
        }
        if(i == nw) fprintf(stderr, "[E::%s] ERROR at %s:%d\n", __func__, __FILE__, __LINE__);
        kmp = kmp | (uint64_t)(aw[i].ol);


        ///here kmp is ul
        kv_push(uint64_t, array->a, kmp);
        kmp = av[k].ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(av[k].v);
        ///here kmp is ol + v
        kv_push(uint64_t, array->a, kmp);
        ///fprintf(stderr, "#av[%u].v>>1: %u, v&1: %u, ol: %u\n", k, av[k].v>>1, av[k].v&1, av[k].ol);
    }

    
    ma_utg_t* tmp;
    for (i = 0; i < realLen; i++)
    {
        w = array->a.a[i];
        tmp = &(ug->u.a[w>>1]);
        if(tmp->m != 0)
        {
            tmp->circ = tmp->end = tmp->len = tmp->m = tmp->n = tmp->start = 0;
            free(tmp->a);
            tmp->a = NULL;
        }
        asg_seq_del(ug->g, w>>1);
    }
    
    ug->u.a[beg_uid>>1] = result;
    ug->g->seq[beg_uid>>1].del = 0;

    uint32_t oLen = 0;
    for (; i < array->a.n; i += 2)
    {
        v = array->a.a[i]>>32;
        w = (uint32_t)array->a.a[i+1];
        /****************************may have bugs********************************/
        ///may have bug here, if there is an edge between beg_uid and end_uid
        ///if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) continue;
        if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) w = v;
        /****************************may have bugs********************************/

        oLen = array->a.a[i+1]>>32;
        asg_append_edges_to_srt(ug->g, v, ug->u.a[v>>1].len, w, oLen, 0, 0, 0);
        oLen = (uint32_t)array->a.a[i];
        asg_append_edges_to_srt(ug->g, w^1, ug->u.a[w>>1].len, v^1, oLen, 0, 0, 0);
    }
}


void init_Edge_iter(asg_t* g, uint32_t v, asg_arc_t* new_edges, uint32_t new_edges_n, Edge_iter* x)
{
	x->av_i = x->new_edges_i = 0;
	x->g = g;
	x->av = asg_arc_a(g, v);
    x->nv = asg_arc_n(g, v);

	if(new_edges == NULL || new_edges_n == 0)
	{
		x->new_edges = NULL;
		x->new_edges_n = 0;
	}
	else
	{
		x->new_edges = new_edges;
		x->new_edges_n = new_edges_n;
	}
}

int get_arc_t(Edge_iter* x, asg_arc_t* get)
{
    for (; x->av_i < x->nv; x->av_i++)
    {
        if(x->av[x->av_i].del) continue;
        get = &(x->av[x->av_i]);
        x->av_i++;
        return 1;
    }


    for (; x->new_edges_i < x->new_edges_n; x->new_edges_i++)
    {
        if(x->new_edges[x->new_edges_i].del) continue;
        get = &(x->new_edges[x->new_edges_i]);
        x->new_edges_i++;
        return 1;
    }

    return 0;
}

///need to consider circle
void merge_ug_nodes_advance(ma_ug_t *ug, asg_t* read_g, R_to_U* ruIndex, kvec_t_u64_warp* array,
asg_arc_t* new_edges, uint32_t new_edges_n)
{
    uint32_t beg_uid = array->a.a[0];
    uint32_t end_uid = array->a.a[array->a.n-1];
    uint32_t realLen = array->a.n;
    uint32_t new_uid = array->a.a[0]>>1;
    uint64_t kmp;



    if(array->a.n == 0) return;
    uint32_t i, v, w;
    ma_utg_t result;
    memset(&result, 0, sizeof(ma_utg_t));
    v = w = (uint32_t)-1;
    for (i = 0; i < array->a.n; i++)
    {
        w = array->a.a[i];
        collect_ma_utg_ts(ug, v, w, &result);
        v = w;
    }

    

    if(result.n == 0) return;
    uint32_t l;
    result.len = 0;
    Edge_iter iter_x, iter_y;
    asg_arc_t get_x, get_y;
    get_x.v=get_x.ul=get_x.ol=0;
    get_y.v=get_y.ul=get_y.ol=0;
    for (i = 0; i < result.n - 1; i++)
    {
        l = (uint32_t)-1;
        v = (uint64_t)(result.a[i])>>32;
        w = (uint64_t)(result.a[i + 1])>>32;
        init_Edge_iter(read_g, v, new_edges, new_edges_n, &iter_x);
        while(get_arc_t(&iter_x, &get_x))
        {
            if(get_x.v == w)
            {
                l = asg_arc_len(get_x);
                break;
            }
        }

        if(l==(uint32_t)-1) fprintf(stderr ,"******error, i: %u\n", i);

        result.a[i] = v; result.a[i] = result.a[i]<<32; result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;

        set_R_to_U(ruIndex, v>>1, new_uid, 1);
    }
    

    if(i < result.n)
    {
        v = (uint64_t)(result.a[i])>>32;
        l = read_g->seq[v>>1].len;
        result.a[i] = v;
        result.a[i] = result.a[i]<<32;
        result.a[i] = result.a[i] | (uint64_t)(l);
        result.len += l;

        set_R_to_U(ruIndex, v>>1, new_uid, 1);
    }
    //has already set result.a, result.len, result.n, result.m
    result.circ = 0;
    result.start = result.a[0]>>32;
    result.end = (result.a[result.n-1]>>32)^1;

    
    


    ///corresponding to direction 1 of new node
    v = beg_uid^1;
    init_Edge_iter(ug->g, v, new_edges, new_edges_n, &iter_x);
    while(get_arc_t(&iter_x, &get_x))
    {
        
        ///if((av[k].v>>1) == (end_uid>>1)) continue;
        w = get_x.v^1;
        init_Edge_iter(ug->g, w, new_edges, new_edges_n, &iter_y);
        while (get_arc_t(&iter_y, &get_y))
        {
            if(get_y.v==(v^1)) break;
        }

        if((get_y.ul>>32)!=w||get_y.v==(v^1)) fprintf(stderr ,"******error\n");
        
        ///here kmp is ul
        kmp = new_uid<<1; kmp = kmp^1; kmp = kmp << 32; kmp=kmp|(uint64_t)(get_y.ol);
        kv_push(uint64_t, array->a, kmp);


        ///here kmp is ol + v
        kmp = get_x.ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(get_x.v);
        kv_push(uint64_t, array->a, kmp);
    }



    ///corresponding to direction 0 of new node
    v = end_uid;
    init_Edge_iter(ug->g, v, new_edges, new_edges_n, &iter_x);
    while(get_arc_t(&iter_x, &get_x))
    {
        
        ///if((av[k].v>>1) == (end_uid>>1)) continue;
        w = get_x.v^1;
        init_Edge_iter(ug->g, w, new_edges, new_edges_n, &iter_y);
        while (get_arc_t(&iter_y, &get_y))
        {
            if(get_y.v==(v^1)) break;
        }

        if((get_y.ul>>32)!=w||get_y.v==(v^1)) fprintf(stderr ,"******error\n");
        
        ///here kmp is ul
        kmp = new_uid<<1; kmp = kmp << 32; kmp=kmp|(uint64_t)(get_y.ol);
        kv_push(uint64_t, array->a, kmp);


        ///here kmp is ol + v
        kmp = get_x.ol; kmp = kmp<<32; kmp = kmp|(uint64_t)(get_x.v);
        kv_push(uint64_t, array->a, kmp);
    }



    
    ma_utg_t* tmp;
    for (i = 0; i < realLen; i++)
    {
        w = array->a.a[i];
        tmp = &(ug->u.a[w>>1]);
        if(tmp->m != 0)
        {
            tmp->circ = tmp->end = tmp->len = tmp->m = tmp->n = tmp->start = 0;
            free(tmp->a);
            tmp->a = NULL;
        }
        asg_seq_del(ug->g, w>>1);
    }
    
    ug->u.a[beg_uid>>1] = result;
    ug->g->seq[beg_uid>>1].del = 0;

    uint32_t oLen = 0;
    for (; i < array->a.n; i += 2)
    {
        v = array->a.a[i]>>32;
        w = (uint32_t)array->a.a[i+1];
        /****************************may have bugs********************************/
        ///may have bug here, if there is an edge between beg_uid and end_uid
        ///if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) continue;
        if(((w>>1) == (beg_uid>>1)) || ((w>>1) == (end_uid>>1))) w = v;
        /****************************may have bugs********************************/

        oLen = array->a.a[i+1]>>32;
        asg_append_edges_to_srt(ug->g, v, ug->u.a[v>>1].len, w, oLen, 0, 0, 0);
        oLen = (uint32_t)array->a.a[i];
        asg_append_edges_to_srt(ug->g, w^1, ug->u.a[w>>1].len, v^1, oLen, 0, 0, 0);
    }


    ///check if it is a circle
    if(get_real_length(ug->g, new_uid, NULL)!=1||get_real_length(ug->g, new_uid^1, NULL)!=1) return;

    get_real_length(ug->g, new_uid, &w);
    if(w!=new_uid) return;
    
    new_uid = new_uid^1;
    get_real_length(ug->g, new_uid, &w);
    if(w!=new_uid) return;

    l = (uint32_t)-1;
    v = (uint64_t)(result.a[result.n - 1])>>32;
    w = (uint64_t)(result.a[0])>>32; 
    init_Edge_iter(read_g, v, new_edges, new_edges_n, &iter_x);
    while(get_arc_t(&iter_x, &get_x))
    {
        if(get_x.v == w)
        {
            l = asg_arc_len(get_x);
            break;
        }
    }

    if(l==(uint32_t)-1) fprintf(stderr ,"******error, i: %u\n", i);

    result.circ = 1; result.start = result.end = UINT32_MAX;
    result.len = result.len - (uint32_t)(result.a[result.n - 1]); result.len = result.len + l;
    result.a[result.n - 1] = v; result.a[result.n - 1] = result.a[result.n - 1]<<32; 
    result.a[result.n - 1] = result.a[result.n - 1] | (uint64_t)(l);


    ug->u.a[beg_uid>>1] = result;
    ug->g->seq[beg_uid>>1].del = 0;
}

void unroll_simple_case(ma_ug_t *ug, asg_t* read_g)
{
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, rnw, nw, w1, w2, beg, end, i;
    asg_arc_t *aw;
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);

    for (v = 0; v < n_vtx; ++v) 
    {
        ///if (nsg->seq[v>>1].del || nsg->seq[v>>1].c == ALTER_LABLE) continue;
        if (nsg->seq[v>>1].del) continue;
        if(asg_arc_n(nsg, v) < 1 || asg_arc_n(nsg, v^1) < 1) continue;
        if(get_real_length(nsg, v, NULL) != 1 || get_real_length(nsg, v^1, NULL) != 1) continue;
        get_real_length(nsg, v, &w1); get_real_length(nsg, v^1, &w2);
        if((v>>1) == (w1>>1) || (v>>1) == (w2>>1)) continue;

        ///for simple circle
        if(w1 == (w2^1))
        {
            beg = end = 0;
            aw = asg_arc_a(nsg, w1^1);
            nw = asg_arc_n(nsg, w1^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == (v^1)) continue;
                beg = aw[i].v;
            }
            if(rnw != 2) continue;

            aw = asg_arc_a(nsg, w2^1);
            nw = asg_arc_n(nsg, w2^1);
            for (i = 0, rnw = 0; i < nw; i++)
            {
                if(aw[i].del) continue;
                rnw++;
                if(aw[i].v == v) continue;
                end = aw[i].v;
            }
            if(rnw != 2) continue;

            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            if((beg>>1) == (end>>1)) continue;
            ///fprintf(stderr, "\n***v>>1: %u\n", v>>1);
            u_vecs.a.n = 0;
            kv_push(uint64_t, u_vecs.a, beg^1);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, v);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, end);
            merge_ug_nodes(ug, read_g, &u_vecs);
        }
        else if(w1 == w2)
        {
            if(get_real_length(nsg, w1^1, NULL) != 2) continue;
            if(get_real_length(nsg, w1, NULL) != 2) continue;
            end = beg = (uint32_t)-1;

            aw = asg_arc_a(nsg, w1);
            nw = asg_arc_n(nsg, w1);
            for (i = 0, rnw = 0; i < nw && rnw < 2; i++)
            {
                if(aw[i].del) continue;
                if(rnw == 0) beg = aw[i].v;
                if(rnw == 1) end = aw[i].v;
                rnw++;
            }
            if((beg>>1) == (end>>1)) continue;
            if(get_real_length(nsg, beg^1, NULL)!=1) continue;
            if(get_real_length(nsg, end^1, NULL)!=1) continue;
            ///fprintf(stderr, "\n###v>>1: %u\n", v>>1);
            u_vecs.a.n = 0;
            kv_push(uint64_t, u_vecs.a, beg^1);
            kv_push(uint64_t, u_vecs.a, w1^1);
            kv_push(uint64_t, u_vecs.a, v);
            kv_push(uint64_t, u_vecs.a, w1);
            kv_push(uint64_t, u_vecs.a, end);
            merge_ug_nodes(ug, read_g, &u_vecs);
        }
    }

    kv_destroy(u_vecs.a);
}

void unroll_simple_case_advance(ma_ug_t *ug, asg_t* read_g, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    asg_t* nsg = ug->g;
    uint32_t v, n_vtx = nsg->n_seq * 2, rnw, nw, beg, end, i;
    uint32_t v_left, v_right, w_left, w_term, w_right, return_flag, convex, /**is_found,**/ n_reduce = 1;
    long long ll, tmp, max_stop_nodeLen, max_stop_baseLen;
    asg_arc_t *aw;
    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);
    buf_t b_0, b_1;
    memset(&b_0, 0, sizeof(buf_t));
    memset(&b_1, 0, sizeof(buf_t));


    while (n_reduce > 0)
    {
        n_reduce = 0;
        ///break nearly circle, forget why...
        n_reduce += asg_arc_del_simple_circle_untig(NULL, NULL, nsg, 100, 0);

        for (v = 0; v < n_vtx; ++v)   // hamt note: using auxsg; read_g is only used when merging unitigs
        {
            if (nsg->seq[v>>1].del) continue;
            v_left = v;
            if(asg_arc_n(nsg, v_left)<1) continue;
            if(get_real_length(nsg, v_left, NULL)!=1) continue;

            return_flag = get_unitig(nsg, NULL, v_left^1, &v_right, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, NULL);
            if(return_flag == LOOP) continue;
            if(return_flag != MUL_INPUT) continue;
            if(asg_arc_n(nsg, v_right)<1) continue;
            if(get_real_length(nsg, v_right, NULL)!=1) continue;

            get_real_length(nsg, v_left, &w_left);
            get_real_length(nsg, v_right, &w_right);
            // forbid circling
            if((v_left>>1)==(w_left>>1)||(v_left>>1)==(w_right>>1)) continue;
            if((v_right>>1)==(w_left>>1)||(v_right>>1)==(w_right>>1)) continue;

            if(w_left!=w_right)
            {
                return_flag = get_unitig(nsg, NULL, w_left, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, NULL);
                if(return_flag == LOOP) continue;
                if(return_flag != MUL_OUTPUT) continue;
                if(convex != (w_right^1)) continue;

                // hamt note: if we reach here, w_left's single path loops all the way back to w_right^1

                beg = end = (uint32_t)-1;
                aw = asg_arc_a(nsg, w_left^1);
                nw = asg_arc_n(nsg, w_left^1);
                for (i = 0, rnw = 0; i < nw; i++)
                {
                    if(aw[i].del) continue;
                    rnw++;
                    if(aw[i].v == (v_left^1)) continue;
                    beg = aw[i].v;
                }
                if(rnw != 2) continue;

                aw = asg_arc_a(nsg, w_right^1);
                nw = asg_arc_n(nsg, w_right^1);
                for (i = 0, rnw = 0; i < nw; i++)
                {
                    if(aw[i].del) continue;
                    rnw++;
                    if(aw[i].v == (v_right^1)) continue;
                    end = aw[i].v;
                }
                if(rnw != 2) continue;

                if(get_real_length(nsg, beg^1, NULL)!=1) continue;
                if(get_real_length(nsg, end^1, NULL)!=1) continue;
                if((beg>>1) == (end>>1)) continue;

                u_vecs.a.n = 0;
                kv_push(uint64_t, u_vecs.a, beg^1);

                b_0.b.n = 0;
                get_unitig(nsg, NULL, w_left, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }
                b_0.b.n = 0;
                get_unitig(nsg, NULL, v_right^1, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }
                b_0.b.n = 0;
                get_unitig(nsg, NULL, w_left, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }

                kv_push(uint64_t, u_vecs.a, end);
                merge_ug_nodes(ug, read_g, &u_vecs);
                n_reduce++;

                // if(nsg->seq[v_left>>1].c!= ALTER_LABLE)
                // {
                //     fprintf(stderr, "Case1: v_left: %u, v_right: %u, w_left: %u, w_right: %u\n",
                //     v_left>>1, v_right>>1, w_left>>1, w_right>>1);
                // }
            }
            else ///if(w_left == w_right)
            {
                if(get_real_length(nsg, w_left^1, NULL)!=2) continue;
                return_flag = get_unitig(nsg, NULL, w_left, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, NULL);
                if(return_flag == LOOP) continue;
                if(return_flag != MUL_OUTPUT) continue;
                if(get_real_length(nsg, convex, NULL)!=2) continue;
                beg = end = (uint32_t)-1;

                aw = asg_arc_a(nsg, convex);
                nw = asg_arc_n(nsg, convex);
                for (i = 0, rnw = 0; i < nw; i++)
                {
                    if(aw[i].del) continue;
                    if(rnw == 0) beg = aw[i].v;
                    if(rnw == 1) end = aw[i].v;
                    rnw++;
                }

                if(rnw != 2) continue;
                if((beg>>1) == (end>>1)) continue;
                if(get_real_length(nsg, beg^1, NULL)!=1) continue;
                if(get_real_length(nsg, end^1, NULL)!=1) continue;

                if(check_different_haps(nsg, ug, read_g, beg, end, reverse_sources, &b_0, &b_1,
                    ruIndex, 2, 1) == PLOID)
                {
                    continue;
                }

                w_term = convex^1;
                u_vecs.a.n = 0;
                kv_push(uint64_t, u_vecs.a, beg^1);

                b_0.b.n = 0;
                get_unitig(nsg, NULL, w_term, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }

                b_0.b.n = 0;
                get_unitig(nsg, NULL, v_left^1, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }

                b_0.b.n = 0;
                get_unitig(nsg, NULL, w_left, &convex, &ll, &tmp, &max_stop_nodeLen, 
                            &max_stop_baseLen, 1, &b_0);
                for (i = 0; i < b_0.b.n; i++)
                {
                    kv_push(uint64_t, u_vecs.a, b_0.b.a[i]);
                }

                kv_push(uint64_t, u_vecs.a, end);
                merge_ug_nodes(ug, read_g, &u_vecs);
                n_reduce++;

                // if(nsg->seq[v_left>>1].c!= ALTER_LABLE)
                // {
                //     fprintf(stderr, "Case2: v_left: %u, v_right: %u, w_left: %u, w_right: %u\n",
                //     v_left>>1, v_right>>1, w_left>>1, w_right>>1);
                // }
            }
        }

        
    }

    free(b_0.b.a);
    free(b_1.b.a);
    kv_destroy(u_vecs.a);
}



void adjust_utg_advance(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    double startTime = Get_T();

    asg_t* nsg = ug->g;
    unroll_simple_case_advance(ug, sg, reverse_sources, ruIndex);
    ///debug_utg_graph(ug, sg, 0, 0);
    drop_semi_circle(ug, ug->g, sg, reverse_sources, ruIndex);
    asg_cleanup(nsg);
    asg_symm(nsg);
    ///debug_utg_graph(ug, sg, 0, 0);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n", __func__, Get_T()-startTime);
    }
}





asg_t *copy_graph(asg_t* src, int round)
{
    asg_t *rg;
	///just calloc
	rg = asg_init();
    uint64_t i, k;
    for (i = 0; i < src->n_seq; ++i) 
    {
        ///if a read has been deleted, should we still add them?
		asg_seq_set(rg, i, src->seq[i].len, src->seq[i].del);
        rg->seq[i].c = src->seq[i].c;
	}
    asg_cleanup(rg);

    
    uint32_t v, n_vtx = src->n_seq * 2, totalL = 0;;


    for (k = 0; k < (uint32_t)round; k++)
    {
        for (i = k; i < src->n_arc; i = i + round)
        {
            if(totalL%50000==0) fprintf(stderr, "totalL: %u, n_arc: %u\n", totalL, src->n_arc);
            totalL++;
            if(src->arc[i].del) continue;
            asg_append_edges_to_srt(rg, src->arc[i].ul>>32, 
            (uint32_t)(src->arc[i].ul) + src->arc[i].ol,
            src->arc[i].v, src->arc[i].ol, src->arc[i].strong, 
            src->arc[i].el, src->arc[i].no_l_indel);
        }
    }
    
    
    totalL = 0;
    for (v = 0; v < n_vtx; ++v)
    {
        if (src->seq[v>>1].del) continue;
        uint32_t nv = asg_arc_n(src, v);
        asg_arc_t *av = asg_arc_a(src, v);
        for (i = 0; i < nv; i++)
        {
            if(totalL%50000==0) fprintf(stderr, "-totalL: %u, n_arc: %u\n", totalL, src->n_arc);
            totalL++;
            if(av[i].del) continue;
            asg_append_edges_to_srt(rg, av[i].ul>>32, (uint32_t)(av[i].ul) + av[i].ol,
            av[i].v, av[i].ol, av[i].strong, av[i].el, av[i].no_l_indel);
        }
    }
    

    

    /**
    for (v = 0; v < n_vtx; ++v)
    {
        if(src->idx[v] != rg->idx[v])
        {
            fprintf(stderr, "*****v: %u, src_i: %u, srcLen: %u, rg_i: %u, rgLen: %u\n", 
            v, src->idx[v]>>32, (uint32_t)src->idx[v], 
            rg->idx[v]>>32, (uint32_t)rg->idx[v]); 
        }
    }
    **/


    for (v = 0; v < n_vtx; ++v)
    {
        
        if(src->seq[v>>1].del != rg->seq[v>>1].del)
        {
            fprintf(stderr, "[E::%s]ERROR1\n", __func__);
        }

        uint32_t nv = asg_arc_n(src, v);
        asg_arc_t *av = asg_arc_a(src, v);
        ///if(nv != rnv)
        if(get_real_length(src, v, NULL) != get_real_length(rg, v, NULL))
        {
            fprintf(stderr, "[E::%s]ERROR2\n", __func__);
        }

        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            asg_arc_t forward = av[i], backward;
            memset(&backward, 0, sizeof(backward));
            if(get_arc(rg, (forward.ul>>32), forward.v, &backward)!=1)
            {
                fprintf(stderr, "[E::%s]ERROR3\n", __func__);
            }

            if(forward.del != backward.del || forward.el != backward.el || 
               forward.no_l_indel != backward.no_l_indel || forward.ol != backward.ol ||
               forward.strong != backward.strong || forward.ul != backward.ul || forward.v != backward.v)
            {
                fprintf(stderr, "[E::%s]ERROR4\n", __func__);
                fprintf(stderr, "forward, del: %u, el: %u, no_l_indel: %u, ol: %u, strong: %u, ul: %u, v: %u\n",
                (uint32_t)forward.del, (uint32_t)forward.el, (uint32_t)forward.no_l_indel, 
                (uint32_t)forward.ol, (uint32_t)forward.strong,
                (uint32_t)forward.ul, (uint32_t)forward.v);
                fprintf(stderr, "backward, del: %u, el: %u, no_l_indel: %u, ol: %u, strong: %u, ul: %u, v: %u\n",
                (uint32_t)backward.del, (uint32_t)backward.el, (uint32_t)backward.no_l_indel, 
                (uint32_t)backward.ol, (uint32_t)backward.strong,
                (uint32_t)backward.ul, (uint32_t)backward.v);
            }
        }

        
        nv = asg_arc_n(rg, v);
        av = asg_arc_a(rg, v);
        for (i = 0; i < nv; i++)
        {
            if(av[i].del) continue;
            asg_arc_t forward = av[i], backward;
            memset(&backward, 0, sizeof(backward));
            if(get_arc(src, (forward.ul>>32), forward.v, &backward)!=1)
            {
                fprintf(stderr, "[E::%s]ERROR5\n", __func__);
            }

            if(forward.del != backward.del || forward.el != backward.el || 
               forward.no_l_indel != backward.no_l_indel || forward.ol != backward.ol ||
               forward.strong != backward.strong || forward.ul != backward.ul || forward.v != backward.v)
            {
                fprintf(stderr, "[E::%s]ERROR6\n", __func__);
            }
        }
        
        
        ///get_arc(g, forward.v^1, (forward.ul>>32)^1, &backward) 
        ///fprintf(stderr, "+v: %u\n", v);
        
    }

    
    if(src->seq_vis)
    {
        rg->seq_vis = (uint8_t*)malloc(src->n_seq*2*sizeof(uint8_t));
        memcpy(rg->seq_vis, src->seq_vis, src->n_seq*2*sizeof(uint8_t));
    }
    

    fprintf(stderr, "end_copy\n");
    return rg;
}

int load_coverage_cut(ma_sub_t** coverage_cut, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    int f_flag = 0;
    uint64_t n_read;
    f_flag += fread(&n_read, sizeof(n_read), 1, fp);
    (*coverage_cut) = (ma_sub_t*)malloc(sizeof(ma_sub_t)*n_read);

    uint64_t i = 0, tmp;
    for (i = 0; i < n_read; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].c = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].del = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].e = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*coverage_cut)[i].s = tmp;
    }
    free(index_name);    
    fflush(fp);
    fclose(fp);
    return 1;
}


int write_coverage_cut(ma_sub_t* coverage_cut, char* read_file_name, uint64_t n_read)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    fwrite(&n_read, sizeof(n_read), 1, fp);
    uint64_t i = 0, tmp;
    for (i = 0; i < n_read; i++)
    {
        tmp = coverage_cut[i].c;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].e;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = coverage_cut[i].s;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int write_ruIndex(R_to_U* ruIndex, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    fwrite(&ruIndex->len, sizeof(ruIndex->len), 1, fp);
    fwrite(ruIndex->index, sizeof(ruIndex->index[0]), ruIndex->len, fp);
    fwrite(R_INF.trio_flag, sizeof(R_INF.trio_flag[0]), ruIndex->len, fp);

    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int load_ruIndex(R_to_U* ruIndex, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    int f_flag = 0;
    f_flag += fread(&(ruIndex)->len, sizeof((ruIndex)->len), 1, fp);
    (ruIndex)->index = (uint32_t*)malloc(sizeof(uint32_t)*(ruIndex)->len);
    f_flag += fread((ruIndex)->index, sizeof((ruIndex)->index[0]), (ruIndex)->len, fp);

    R_INF.trio_flag = (uint8_t*)malloc(sizeof(uint8_t)*(ruIndex)->len);
    f_flag += fread(R_INF.trio_flag, sizeof(R_INF.trio_flag[0]), (ruIndex)->len, fp);

    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}

int write_asg_t(asg_t *sg, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "w");
    uint32_t tmp, i;
    
    tmp = sg->n_arc;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    // tmp = sg->m_arc;
    fwrite(&tmp, sizeof(tmp), 1, fp);

    tmp = sg->is_srt;
    fwrite(&tmp, sizeof(tmp), 1, fp);


    tmp = sg->n_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    // tmp = sg->m_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    
    tmp = sg->is_symm;
    fwrite(&tmp, sizeof(tmp), 1, fp);
    tmp = sg->r_seq;
    fwrite(&tmp, sizeof(tmp), 1, fp);


    uint32_t Len;

    Len = sg->n_seq*2;
    fwrite(sg->seq_vis, sizeof(sg->seq_vis[0]), Len, fp);

    Len = sg->n_seq*2;
    fwrite(sg->idx, sizeof(sg->idx[0]), Len, fp);


    for (i = 0; i < sg->n_arc; i++)
    {
        tmp = sg->arc[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].el;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].no_l_indel;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].ol;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->arc[i].strong;
        fwrite(&tmp, sizeof(tmp), 1, fp);

        uint64_t tmp_64;
        tmp_64 = sg->arc[i].ul;
        fwrite(&tmp_64, sizeof(tmp_64), 1, fp);

        tmp = sg->arc[i].v;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }


    for (i = 0; i < sg->n_seq; i++)
    {
        tmp = sg->seq[i].c;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->seq[i].del;
        fwrite(&tmp, sizeof(tmp), 1, fp);
        tmp = sg->seq[i].len;
        fwrite(&tmp, sizeof(tmp), 1, fp);
    }
    
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}


int load_asg_t(asg_t **sg, char* read_file_name)
{
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.bin", read_file_name);
    FILE* fp = fopen(index_name, "r");
    if(!fp)
    {
        return 0;
    }
    uint32_t tmp, i;

    (*sg) = (asg_t*)calloc(1, sizeof(asg_t));
    int f_flag = 0;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->n_arc = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->m_arc = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->is_srt = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->n_seq = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->m_seq = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->is_symm = tmp;
    f_flag += fread(&tmp, sizeof(tmp), 1, fp);
    (*sg)->r_seq = tmp;

    uint32_t Len;

    Len = (*sg)->n_seq*2;
    (*sg)->seq_vis = (uint8_t*)malloc(sizeof(uint8_t)*Len);
    f_flag += fread((*sg)->seq_vis, sizeof((*sg)->seq_vis[0]), Len, fp);



    Len = (*sg)->n_seq*2;
    (*sg)->idx = (uint64_t*)malloc(sizeof(uint64_t)*Len);
    f_flag += fread((*sg)->idx, sizeof((*sg)->idx[0]), Len, fp);
    
    



    (*sg)->arc = (asg_arc_t*)malloc(sizeof(asg_arc_t)*(*sg)->m_arc);
    (*sg)->seq = (asg_seq_t*)malloc(sizeof(asg_seq_t)*(*sg)->m_seq);
    

    for (i = 0; i < (*sg)->n_arc; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].del = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].el = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].no_l_indel = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].ol = tmp;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].strong = tmp;

        uint64_t tmp_64;
        f_flag += fread(&tmp_64, sizeof(tmp_64), 1, fp);
        (*sg)->arc[i].ul = tmp_64;

        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->arc[i].v = tmp;
    }


    for (i = 0; i < (*sg)->n_seq; i++)
    {
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].c = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].del = tmp;
        f_flag += fread(&tmp, sizeof(tmp), 1, fp);
        (*sg)->seq[i].len = tmp;
    }
    
    free(index_name);    
    fflush(fp);
    fclose(fp);

    return 1;
}


int write_debug_graph(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{

    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);

    ////write_All_reads(&R_INF, gfa_name);
    sprintf(gfa_name, "%s.all.debug.source", output_file_name);
    write_ma_hit_ts(sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.debug.reverse", output_file_name);
    write_ma_hit_ts(reverse_sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.debug.coverage_cut", output_file_name);
    write_coverage_cut(coverage_cut, gfa_name, n_read);
    sprintf(gfa_name, "%s.all.debug.ruIndex", output_file_name);
    write_ruIndex(ruIndex, gfa_name);
    sprintf(gfa_name, "%s.all.debug.asg_t", output_file_name);
    write_asg_t(sg, gfa_name);
    free(gfa_name);
    fprintf(stderr, "debug_graph has been written.\n");
    return 1;
}
int hamt_write_debug_graph(asg_t *sg, ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, 
char* output_file_name, long long n_read, ma_hit_t_alloc* reverse_sources, R_to_U* ruIndex)
{
    // monkey patch of a debug utility: set probing point in the graph cleaning flow and 
    //                               write/read graph bins from there (instead of starting from 
    //                               raw unitig graph every time.)

    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);

    ////write_All_reads(&R_INF, gfa_name);
    sprintf(gfa_name, "%s.all.probe.source", output_file_name);
    write_ma_hit_ts(sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.probe.reverse", output_file_name);
    write_ma_hit_ts(reverse_sources, R_INF.total_reads, gfa_name);
    sprintf(gfa_name, "%s.all.probe.coverage_cut", output_file_name);
    write_coverage_cut(coverage_cut, gfa_name, n_read);
    sprintf(gfa_name, "%s.all.probe.ruIndex", output_file_name);
    write_ruIndex(ruIndex, gfa_name);
    sprintf(gfa_name, "%s.all.probe.asg_t", output_file_name);
    write_asg_t(sg, gfa_name);
    free(gfa_name);
    fprintf(stderr, "Wrote probing dbg-gfa.\n");
    return 1;
}
int hamt_load_debug_graph(asg_t** sg, ma_hit_t_alloc** sources, ma_sub_t** coverage_cut, 
                          char* output_file_name, 
                          ma_hit_t_alloc** reverse_sources, 
                          R_to_U* ruIndex, long long total_reads)
{
    FILE* fp = NULL;
    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);

    sprintf(gfa_name, "%s.all.probe.source.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {fprintf(stderr, "[W::%s] do not have %s\n", __func__, gfa_name);free(gfa_name); return 0;}
    fclose(fp);

    sprintf(gfa_name, "%s.all.probe.reverse.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {fprintf(stderr, "[W::%s] do not have %s\n", __func__, gfa_name);free(gfa_name); return 0;}
    fclose(fp);

    sprintf(gfa_name, "%s.all.probe.coverage_cut.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {fprintf(stderr, "[W::%s] do not have %s\n", __func__, gfa_name);free(gfa_name); return 0;}
    fclose(fp);

    sprintf(gfa_name, "%s.all.probe.ruIndex.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {fprintf(stderr, "[W::%s] do not have %s\n", __func__, gfa_name);free(gfa_name); return 0;}
    fclose(fp);

    sprintf(gfa_name, "%s.all.probe.asg_t.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {fprintf(stderr, "[W::%s] do not have %s\n", __func__, gfa_name);free(gfa_name); return 0;}
    fclose(fp);

    if((*sources)!=NULL) destory_ma_hit_t_alloc((*sources), total_reads);
    if((*reverse_sources)!=NULL) destory_ma_hit_t_alloc((*reverse_sources), total_reads);
    if((*coverage_cut)!=NULL) free((*coverage_cut));
    if((ruIndex)!=NULL) destory_R_to_U((ruIndex));
    if((*sg)!=NULL) asg_destroy(*sg);

    sprintf(gfa_name, "%s.all.probe.source", output_file_name);
    if (!load_ma_hit_ts(sources, gfa_name)){
        fprintf(stderr, "[W::%s] have %s but loading failed\n", __func__, gfa_name);
        free(gfa_name); return 0;
    }

    sprintf(gfa_name, "%s.all.probe.reverse", output_file_name);
    if (!load_ma_hit_ts(reverse_sources, gfa_name)){
        fprintf(stderr, "[W::%s] have %s but loading failed\n", __func__, gfa_name);
        free(gfa_name); return 0;
    }
    
    sprintf(gfa_name, "%s.all.probe.coverage_cut", output_file_name);
    if (!load_coverage_cut(coverage_cut, gfa_name)){
        fprintf(stderr, "[W::%s] have %s but loading failed\n", __func__, gfa_name);
        free(gfa_name); return 0;
    }

    sprintf(gfa_name, "%s.all.probe.ruIndex", output_file_name);
    if (!load_ruIndex(ruIndex, gfa_name)){
        fprintf(stderr, "[W::%s] have %s but loading failed\n", __func__, gfa_name);
        free(gfa_name); return 0;
    }

    sprintf(gfa_name, "%s.all.probe.asg_t", output_file_name);
    if (!load_asg_t(sg, gfa_name)){
        fprintf(stderr, "[W::%s] have %s but loading failed\n", __func__, gfa_name);
        free(gfa_name); return 0;
    }

    R_INF.paf = (*sources); R_INF.reverse_paf = (*reverse_sources);

    free(gfa_name);
    return 2;
}



int load_debug_graph(asg_t** sg, ma_hit_t_alloc** sources, ma_sub_t** coverage_cut, 
char* output_file_name, ma_hit_t_alloc** reverse_sources, R_to_U* ruIndex, long long total_reads)
{
    FILE* fp = NULL;
    char* gfa_name = (char*)malloc(strlen(output_file_name)+55);
    sprintf(gfa_name, "%s.all.debug.source.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {free(gfa_name); return 0;}
    fclose(fp);
    sprintf(gfa_name, "%s.all.debug.reverse.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {free(gfa_name); return 0;}
    fclose(fp);
    sprintf(gfa_name, "%s.all.debug.coverage_cut.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {free(gfa_name); return 0;}
    fclose(fp);
    sprintf(gfa_name, "%s.all.debug.ruIndex.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {free(gfa_name); return 0;}
    fclose(fp);
    sprintf(gfa_name, "%s.all.debug.asg_t.bin", output_file_name);
    fp = fopen(gfa_name, "r"); if(!fp) {free(gfa_name); return 0;}
    fclose(fp);


    if((*sources)!=NULL)
    {
        destory_ma_hit_t_alloc((*sources), total_reads);
    }

    if((*reverse_sources)!=NULL)
    {
        destory_ma_hit_t_alloc((*reverse_sources), total_reads);
    }

    if((*coverage_cut)!=NULL)
    {
        free((*coverage_cut));
    }

    if((ruIndex)!=NULL)
    {
        destory_R_to_U((ruIndex));
    }

    if((*sg)!=NULL)
    {
        asg_destroy(*sg);
    }

    

    sprintf(gfa_name, "%s.all.debug.source", output_file_name);
    if(!load_ma_hit_ts(sources, gfa_name))
    {
        free(gfa_name);
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.reverse", output_file_name);
    if(!load_ma_hit_ts(reverse_sources, gfa_name))
    {
        free(gfa_name);
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.coverage_cut", output_file_name);
    if(!load_coverage_cut(coverage_cut, gfa_name))
    {
        free(gfa_name);
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.ruIndex", output_file_name);
    if(!load_ruIndex(ruIndex, gfa_name))
    {
        free(gfa_name);
        return 0;
    }

    sprintf(gfa_name, "%s.all.debug.asg_t", output_file_name);
    if(!load_asg_t(sg, gfa_name))
    {
        free(gfa_name);
        return 0;
    }

    R_INF.paf = (*sources); R_INF.reverse_paf = (*reverse_sources);

    free(gfa_name);
    return 1;
}




void print_utg_hap(ma_ug_t *ug, asg_t* read_g, uint32_t uid, ma_hit_t_alloc* reverse_sources, 
R_to_U* ruIndex)
{
    
    uint32_t k, rId, qn, i, is_Unitig;
    ma_utg_t* u = &(ug->u.a[uid]);
    fprintf(stderr, "\nuid: %u, u->n: %u\n", uid, u->n);
    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        fprintf(stderr, "self[%u]: %u\n", k, rId);
    }

    for (k = 0; k < u->n; k++)
    {
        qn = u->a[k]>>33;
        fprintf(stderr, "****self[%u]: %u, Len: %u\n", k, qn, reverse_sources[qn].length);
        for (i = 0; i < reverse_sources[qn].length; i++)
        {
            rId = Get_tn(reverse_sources[qn].buffer[i]);
            if(read_g->seq[rId].del == 1)
            {
                ///get the id of read that contains it 
                get_R_to_U(ruIndex, rId, &rId, &is_Unitig);
                if(rId == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[rId].del == 1) continue;
            }
            fprintf(stderr, "hap[%u]: %u\n", i, rId);
        }
    }

}


void append_utg(ma_ug_t* ptg, ma_ug_t* atg)
{
    uint64_t num_nodes = 0;
    asg_t* nsg = atg->g;
    uint32_t v, n_vtx = nsg->n_seq;
    ma_utg_t *p;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del || atg->u.a[v].m == 0) continue;
        num_nodes++;
    }

    if(num_nodes == 0) return;

    // prepare memory
    ptg->u.n = ptg->u.n + num_nodes;
    if(ptg->u.n > ptg->u.m)
    {
        ptg->u.m = ptg->u.n;
        ptg->u.a = (ma_utg_t*)realloc(ptg->u.a, ptg->u.m*sizeof(ma_utg_t));     
    }
    ptg->u.n = ptg->u.n - num_nodes;

    // migrate
    for (v = 0; v < atg->g->n_seq; ++v) 
    {
        if(atg->g->seq[v].del || atg->u.a[v].m == 0) continue;

        p = &(ptg->u.a[ptg->u.n]);
        p->len = atg->u.a[v].len;
        p->circ = atg->u.a[v].circ;
        p->start = atg->u.a[v].start;
        p->end = atg->u.a[v].end;
        p->m = atg->u.a[v].m; atg->u.a[v].m = 0;
        p->n = atg->u.a[v].n; atg->u.a[v].n = 0;
        p->a = atg->u.a[v].a; atg->u.a[v].a = 0;
        p->s = atg->u.a[v].s; atg->u.a[v].s = 0;
        asg_seq_set(ptg->g, ptg->u.n, p->len, 0);
        ptg->u.n++;
    }

    if(ptg->g->idx != 0) free(ptg->g->idx); 
    ptg->g->idx = 0;
    asg_cleanup(ptg->g);
}


void print_utg_coverage(ma_ug_t *ug, ma_sub_t* coverage_cut, uint32_t v, ma_hit_t_alloc* sources)
{
    asg_t* nsg = ug->g;
    uint32_t rId, k, j;
    ma_utg_t* u = NULL;
    ma_hit_t *h;

    if(nsg->seq[v].del) return;
    u = &(ug->u.a[v]);
    if(u->m == 0) return;
    long long R_bases = 0, C_bases = 0;
    long long U_R_bases = 0, U_C_bases = 0;
    for (k = 0; k < u->n; k++)
    {
        rId = u->a[k]>>33;
        C_bases = 0;
        R_bases = coverage_cut[rId].e - coverage_cut[rId].s;
        for (j = 0; j < (uint64_t)(sources[rId].length); j++)
        {
            h = &(sources[rId].buffer[j]);
            if(h->el != 1) continue;
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
        U_R_bases += R_bases;
        U_C_bases += C_bases;
        C_bases = C_bases/R_bases;

        fprintf(stderr, "%.*s\t%lld\n", (int)Get_NAME_LENGTH(R_INF, rId), Get_NAME(R_INF, rId), C_bases);
    }

    fprintf(stderr, "v: %u, coverage: %lld\n\n", v, U_C_bases/U_R_bases);
}

void recover_utg_by_coverage(ma_ug_t **ptg, asg_t* read_g, ma_sub_t* coverage_cut, 
ma_hit_t_alloc* sources, R_to_U* ruIndex)
{
    if(asm_opt.recover_atg_cov_min == -1) return;
    if(asm_opt.recover_atg_cov_max == -1) return;
    if(asm_opt.recover_atg_cov_min > asm_opt.recover_atg_cov_max) return;
    ma_ug_t *atg = NULL;
    atg = ma_ug_gen_primary(read_g, ALTER_LABLE);
    asg_t* nsg = atg->g;
    uint32_t v, n_vtx = nsg->n_seq, k, j, rId, available_reads = 0, keep_atg = 0, tn, is_Unitig;
    ma_utg_t* u = NULL;
    ma_hit_t *h;

    
    ///print_untig_by_read(atg, "SRR11606870.634978", -1, NULL, NULL, "debug");
    long long R_bases = 0, C_bases = 0, C_bases_primary = 0, C_bases_alter = 0;
    long long total_C_bases = 0, total_C_bases_alter = 0;
    for (v = 0; v < n_vtx; ++v)  // hamt note: using atg's auxsg
    {
        if(nsg->seq[v].del) continue;
        u = &(atg->u.a[v]);
        if(u->m == 0) continue;
        total_C_bases = total_C_bases_alter = available_reads = 0;
        
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            C_bases = C_bases_primary = C_bases_alter = 0;
            R_bases = coverage_cut[rId].e - coverage_cut[rId].s;
            for (j = 0; j < (uint64_t)(sources[rId].length); j++)
            {
                h = &(sources[rId].buffer[j]);
                if(h->el != 1) continue;
                tn = Get_tn((*h));
                if(read_g->seq[tn].del == 1)
                {
                    ///get the id of read that contains it 
                    get_R_to_U(ruIndex, tn, &tn, &is_Unitig);
                    if(tn == (uint32_t)-1 || is_Unitig == 1 || read_g->seq[tn].del == 1) continue;
                }
                if(read_g->seq[tn].del == 1) continue;
                if(read_g->seq[tn].c == ALTER_LABLE)
                {
                    C_bases_alter += Get_qe((*h)) - Get_qs((*h));
                }
                else
                {
                    C_bases_primary += Get_qe((*h)) - Get_qs((*h));
                }
            }

            C_bases = C_bases_primary + C_bases_alter;
            total_C_bases += C_bases; 

            C_bases = C_bases/R_bases;
            if(C_bases >= asm_opt.recover_atg_cov_min && C_bases <= asm_opt.recover_atg_cov_max)
            {
                if(C_bases_alter >= (C_bases_primary + C_bases_alter) * ALTER_COV_THRES) available_reads++;
                total_C_bases_alter += C_bases_alter;
            }
        }


        if(((available_reads < (u->n * 0.8)) && (total_C_bases_alter < (total_C_bases * 0.8))) 
           || available_reads == 0)
        {
            asg_seq_del(nsg, v);

            if(u->m!=0)
            {
                u->m = u->n = 0;
                free(u->a);
                u->a = NULL;
            }
        }
        else
        {
            ///print_utg_coverage(atg, coverage_cut, v, sources);
            ///fprintf(stderr, "rId: %u\n", rId);
            ///fprintf(stderr, "%.*s\t%lld\n", (int)Get_NAME_LENGTH(R_INF, rId), Get_NAME(R_INF, rId), C_bases);
            keep_atg++;
        }
    }

    if(keep_atg > 0)
    {
        asg_cleanup(nsg);
        asg_symm(nsg);
        append_utg(*ptg, atg);

        n_vtx = read_g->n_seq;
        for (v = 0; v < n_vtx; v++)
        {
            read_g->seq[v].c = ALTER_LABLE;
        }
        // renew sg's c labels with contig graph (ptg)
        nsg = (*ptg)->g;
        n_vtx = nsg->n_seq;
        for (v = 0; v < n_vtx; ++v) 
        {
            if(nsg->seq[v].del) continue;
            u = &((*ptg)->u.a[v]);
            if(u->m == 0) continue;
            for (k = 0; k < u->n; k++)
            {
                rId = u->a[k]>>33;
                read_g->seq[rId].c = nsg->seq[v].c;
            }
        }
        // delete reads that's not in the contig graph
        n_vtx = read_g->n_seq;
        for (v = 0; v < n_vtx; v++)
        {
            if(read_g->seq[v].c == ALTER_LABLE)
            {
                asg_seq_drop(read_g, v);
            }
        }
    }

    ///fprintf(stderr, "keep_atg: %u\n", keep_atg);
    ma_ug_destroy(atg);
}




void adjust_utg_by_primary(ma_ug_t **ug, asg_t* read_g, float drop_rate,
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, ma_sub_t* coverage_cut, 
long long bubble_dist, long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp,
kvec_asg_arc_t_warp* new_rtg_edges)
{
    asg_t* nsg = (*ug)->g;
    uint32_t v, n_vtx = nsg->n_seq, k, rId, just_contain;
    ma_utg_t* u = NULL;
    
    ///print_utg_coverage(*ug, coverage_cut, 440, sources);
    ///exit(0);
    /**
    kvec_t_u32_warp new_rtg_nodes;
    kv_init(new_rtg_nodes.a);
    **/

    drop_semi_circle((*ug), nsg, read_g, reverse_sources, ruIndex);
    asg_cleanup(nsg);
    adjust_utg_advance(read_g, (*ug), reverse_sources, ruIndex); // unroll_simple_case_advance and drop_semi_circle
     
    nsg = (*ug)->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen((*ug)->u, v) = (*ug)->u.a[v].n;
    }

    // bubble, tip and triangle treatments
    clean_primary_untig_graph(*ug, read_g, reverse_sources, bubble_dist, tipsLen, 
    tip_drop_ratio, stops_threshold, ruIndex, NULL, NULL, 0, 0, 0, 
    chimeric_rate, 0, 0, drop_ratio);

    delete_useless_nodes(ug);  // completely remove seuqneces that: nsg->seq[v].c == ALTER_LABLE
    renew_utg(ug, read_g, new_rtg_edges);


    #if 0

    if(asm_opt.purge_level_primary > 0)
    {
        just_contain = 0;
        if(asm_opt.purge_level_primary == 1) just_contain = 1;

        purge_dups(*ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, new_rtg_edges, 
        asm_opt.purge_simi_rate, asm_opt.purge_overlap_len, max_hang, min_ovlp, bubble_dist, 
        drop_ratio, just_contain, 0);
        delete_useless_nodes(ug);
        renew_utg(ug, read_g, new_rtg_edges);
    }


    if (!(asm_opt.flag & HA_F_BAN_POST_JOIN))
    {
        rescue_missing_overlaps_aggressive(*ug, read_g, sources, coverage_cut, ruIndex, max_hang,
        min_ovlp, 0, 0, 1, NULL);
        renew_utg(ug, read_g, new_rtg_edges);

        rescue_contained_reads_aggressive(*ug, read_g, sources, coverage_cut, ruIndex, max_hang, 
        min_ovlp, 0, 10, 0, 1, NULL, NULL);
        renew_utg(ug, read_g, new_rtg_edges);

        if(asm_opt.purge_level_primary > 0)
        {
            just_contain = 0;
            if(asm_opt.purge_level_primary == 1) just_contain = 1;

            purge_dups(*ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, new_rtg_edges, 
            asm_opt.purge_simi_rate, asm_opt.purge_overlap_len, max_hang, min_ovlp, bubble_dist, 
            drop_ratio, just_contain, 0);
            delete_useless_nodes(ug);
            renew_utg(ug, read_g, new_rtg_edges);
        }
    }

    if(asm_opt.purge_level_primary == 0)
    {
        purge_dups(*ug, read_g, coverage_cut, sources, reverse_sources, ruIndex, new_rtg_edges, 
        asm_opt.purge_simi_rate, asm_opt.purge_overlap_len, max_hang, min_ovlp, bubble_dist, 
        drop_ratio, 0, 1);
    }

    n_vtx = read_g->n_seq;
    for (v = 0; v < n_vtx; v++)
    {
        read_g->seq[v].c = ALTER_LABLE;
    }
    
    
    nsg = (*ug)->g;
    n_vtx = nsg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        if(nsg->seq[v].c == ALTER_LABLE) continue;
        u = &((*ug)->u.a[v]);
        if(u->m == 0) continue;
        for (k = 0; k < u->n; k++)
        {
            rId = u->a[k]>>33;
            read_g->seq[rId].c = nsg->seq[v].c;
        }
    }

    n_vtx = read_g->n_seq;
    for (v = 0; v < n_vtx; v++)
    {
        if(read_g->seq[v].c == ALTER_LABLE)
        {
            asg_seq_drop(read_g, v);
        }
    }

    if(asm_opt.recover_atg_cov_min == -1024)
    {
        asm_opt.recover_atg_cov_max = asm_opt.hom_global_coverage/HOM_PEAK_RATE;
        asm_opt.recover_atg_cov_min = asm_opt.recover_atg_cov_max * 0.8;
        ///asm_opt.recover_atg_cov_max = asm_opt.recover_atg_cov_max * 1.2;
        asm_opt.recover_atg_cov_max = INT32_MAX;
    }

    if(asm_opt.recover_atg_cov_max != INT32_MAX)
    {
        fprintf(stderr, "[M::%s] primary contig coverage range: [%d, %d]\n", 
        __func__, asm_opt.recover_atg_cov_min, asm_opt.recover_atg_cov_max);
    }
    else
    {
        fprintf(stderr, "[M::%s] primary contig coverage range: [%d, infinity]\n", 
        __func__, asm_opt.recover_atg_cov_min);
    }

    recover_utg_by_coverage(ug, read_g, coverage_cut, sources, ruIndex);
    #endif
    /**
    kv_destroy(new_rtg_nodes.a);
    **/
}


void output_contig_graph_primary_pre(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long bubble_dist, 
long long tipsLen, R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);

    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);

    asg_t* nsg = ug->g;
    uint32_t n_vtx = nsg->n_seq, v;
    for (v = 0; v < n_vtx; ++v) 
    {
        if(nsg->seq[v].del) continue;
        nsg->seq[v].c = PRIMARY_LABLE;
        EvaluateLen(ug->u, v) = ug->u.a[v].n;
    }
    asg_pop_bubble_primary_trio(ug, bubble_dist, (uint32_t)-1, DROP);
    cut_trio_tip_primary(ug->g, ug, tipsLen, (uint32_t)-1, 0, sg, reverse_sources, ruIndex, 2);
    asg_pop_bubble_primary_trio(ug, bubble_dist, (uint32_t)-1, DROP);
    cut_trio_tip_primary(ug->g, ug, tipsLen, (uint32_t)-1, 0, sg, reverse_sources, ruIndex, 2);
    delete_useless_nodes(&ug);
    renew_utg(&ug, sg, &new_rtg_edges);

    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "Writing processed unitig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.p_utg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);
    
    sprintf(gfa_name, "%s.p_utg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "utg", output_file);
    fclose(output_file);
    if(asm_opt.bed_inconsist_rate != 0)
    {
        sprintf(gfa_name, "%s.p_utg.lowQ.bed", output_file_name);
        output_file = fopen(gfa_name, "w");
        ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
        max_hang, min_ovlp, asm_opt.bed_inconsist_rate, "utg", output_file);
        fclose(output_file);
    }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}

void output_contig_graph_primary(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name, 
ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, long long bubble_dist, 
long long tipsLen, float tip_drop_ratio, long long stops_threshold, 
R_to_U* ruIndex, float chimeric_rate, float drop_ratio, int max_hang, int min_ovlp)
{
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, PRIMARY_LABLE);

    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);


    adjust_utg_by_primary(&ug, sg, TRIO_THRES, sources, reverse_sources, coverage_cut, 
    bubble_dist, tipsLen, tip_drop_ratio, stops_threshold, ruIndex, chimeric_rate, drop_ratio, 
    max_hang, min_ovlp, &new_rtg_edges);

    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);
    
    fprintf(stderr, "Writing primary contig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.p_ctg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "ptg", output_file);
    fclose(output_file);

    sprintf(gfa_name, "%s.p_ctg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "ptg", output_file);
    fclose(output_file);
    if(asm_opt.bed_inconsist_rate != 0)
    {
        sprintf(gfa_name, "%s.p_ctg.lowQ.bed", output_file_name);
        output_file = fopen(gfa_name, "w");
        ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
        max_hang, min_ovlp, asm_opt.bed_inconsist_rate, "ptg", output_file);
        fclose(output_file);
    }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}




void output_contig_graph_alternative(asg_t *sg, ma_sub_t* coverage_cut, char* output_file_name,
ma_hit_t_alloc* sources, R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    kvec_asg_arc_t_warp new_rtg_edges;
    kv_init(new_rtg_edges.a);
    ma_ug_t *ug = NULL;
    ug = ma_ug_gen_primary(sg, ALTER_LABLE);
    ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    fprintf(stderr, "Writing alternate contig GFA to disk... \n");
    char* gfa_name = (char*)malloc(strlen(output_file_name)+35);
    sprintf(gfa_name, "%s.a_ctg.gfa", output_file_name);
    FILE* output_file = fopen(gfa_name, "w");
    ma_ug_print(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "atg", output_file);
    fclose(output_file);
    
    sprintf(gfa_name, "%s.a_ctg.noseq.gfa", output_file_name);
    output_file = fopen(gfa_name, "w");
    ma_ug_print_simple(ug, &R_INF, sg, coverage_cut, sources, ruIndex, "atg", output_file);
    fclose(output_file);
    if(asm_opt.bed_inconsist_rate != 0)
    {
        sprintf(gfa_name, "%s.a_ctg.lowQ.bed", output_file_name);
        output_file = fopen(gfa_name, "w");
        ma_ug_print_bed(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, 
        max_hang, min_ovlp, asm_opt.bed_inconsist_rate, "atg", output_file);
        fclose(output_file);
    }

    free(gfa_name);
    ma_ug_destroy(ug);
    kv_destroy(new_rtg_edges.a);
}

int output_tips(asg_t *g, const All_reads *RNF)
{
	uint32_t v, n_vtx = g->n_seq * 2;
    for (v = 0; v < n_vtx; ++v)
    {
        if (g->seq[v>>1].del) continue;

        if(asg_arc_n(g, v) == 0)
        {
            fprintf(stderr, "%.*s\n",
             (int)Get_NAME_LENGTH((*RNF), v>>1), 
             Get_NAME((*RNF), v>>1));
        }
    }

    return 1;
}

void collect_abnormal_edges(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum)
{
    double startTime = Get_T();
    long long T_edges, T_Single_Dir_Edges_0, T_Single_Dir_Edges_1, T_Conflict_Equal_Edges, T_Conflict_Strong_Edges;
    T_edges = T_Single_Dir_Edges_0 = T_Single_Dir_Edges_1 = T_Conflict_Equal_Edges = T_Conflict_Strong_Edges = 0;
    long long T_Single_Dir_Edges_1_1000 = 0;
    long long related_reads = 0;
    long long related_overlaps = 0;
    long long i, j;
    uint32_t qn, tn;
    int is_equal_f, is_strong_f; 
    int is_equal_b, is_strong_b, is_exist_b; 

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    
    for (i = 0; i < readNum; i++)
    {
        for (j = 0; j < (long long)paf[i].length; j++)
        {
            qn = Get_qn(paf[i].buffer[j]);
            tn = Get_tn(paf[i].buffer[j]);
            T_edges++;

            is_equal_f = paf[i].buffer[j].el;
            is_strong_f = paf[i].buffer[j].ml;
            
            is_exist_b = get_specific_overlap(&(paf[tn]), tn, qn);
            if(is_exist_b == -1)
            {
                is_exist_b = get_specific_overlap(&(rev_paf[tn]), tn, qn);
                if(is_exist_b != -1)
                {
                    T_Single_Dir_Edges_0++;
                    kv_push(uint64_t, edge_vector, qn);
                    kv_push(uint64_t, edge_vector, tn);
                    ///related_overlaps += paf[qn].length + rev_paf[qn].length + paf[tn].length + rev_paf[tn].length;
                    // fprintf(stderr, "%.*s(%d) ---(+)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(-)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(rev_paf[tn].buffer[is_exist_b]) - Get_qs(rev_paf[tn].buffer[is_exist_b]));
                } 
                else
                {
                    T_Single_Dir_Edges_1++;
                    if(Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]) >= 1000)
                    {
                        T_Single_Dir_Edges_1_1000++;
                        // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                        // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                        // is_strong_f,
                        // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                        // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));
                    }
                }

                related_reads = related_reads + 2;
            }
            else
            {
                is_equal_b = paf[tn].buffer[is_exist_b].el;
                is_strong_b = paf[tn].buffer[is_exist_b].ml;

                if(is_equal_f != is_equal_b)
                {
                    T_Conflict_Equal_Edges++;

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // is_equal_f,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // is_equal_b,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(paf[tn].buffer[is_exist_b]) - Get_qs(paf[tn].buffer[is_exist_b]));
                } 
                
                if(is_strong_f != is_strong_b)
                {
                    T_Conflict_Strong_Edges++;
                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n", 
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // is_strong_f,
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // Get_qe(paf[i].buffer[j]) - Get_qs(paf[i].buffer[j]));

                    // fprintf(stderr, "%.*s(%d) ---(%d)--> %.*s(%d), Len: %d\n\n", 
                    // Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn), tn,
                    // is_strong_b,
                    // Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn), qn,
                    // Get_qe(paf[tn].buffer[is_exist_b]) - Get_qs(paf[tn].buffer[is_exist_b]));
                } 

                if(is_equal_f != is_equal_b || is_strong_f != is_strong_b)
                {
                    related_reads++;
                }
            }
        }
    }

    radix_sort_arch64(edge_vector.a, edge_vector.a + edge_vector.n);
    uint64_t pre = (uint64_t)-1;
    long long mn = 0;
    for (i = 0; i < (long long)edge_vector.n; i++)
    {
        if(pre != edge_vector.a[i])
        {
            mn++;
            pre = edge_vector.a[i];
            related_overlaps += paf[pre].length + rev_paf[pre].length;
        }

        if(i>0 && edge_vector.a[i] < edge_vector.a[i-1]) fprintf(stderr, "hehe\n");
    }
    
    

    fprintf(stderr, "****************statistic for abnormal overlaps****************\n");
    fprintf(stderr, "overlaps #: %lld\n", T_edges);
	fprintf(stderr, "one direction overlaps (different phasing)#: %lld\n", T_Single_Dir_Edges_0);
    fprintf(stderr, "one direction overlaps (missing)#: %lld\n", T_Single_Dir_Edges_1);
    fprintf(stderr, "one direction overlaps (missing) >= 1000#: %lld\n", T_Single_Dir_Edges_1_1000);
    fprintf(stderr, "conflict strong/weak overlaps #: %lld\n", T_Conflict_Strong_Edges);
    fprintf(stderr, "conflict exact/inexact overlaps #: %lld\n", T_Conflict_Equal_Edges);
    fprintf(stderr, "related_reads #: %lld/%lld\n", related_reads, mn);
    fprintf(stderr, "related_overlaps #: %lld\n", related_overlaps);
    fprintf(stderr, "****************statistic for abnormal overlaps****************\n");

    fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);

    kv_destroy(edge_vector);
}

///if we don't have this function, we just simply remove all one-direction edges
///by utilizing this function, some one-direction edges can be recovered as two-direction edges
///trio does not influence this function  
void try_rescue_overlaps(ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, long long readNum, long long rescue_threshold)
// called as (sources, reverse_sources, n_read, 4);  where n_read is R_INF.total_reads
{
    double startTime = Get_T();
    long long i, j, revises = 0;
    uint32_t qn, tn, qs, qe;

    kvec_t(uint64_t) edge_vector;
    kv_init(edge_vector);
    kvec_t(uint64_t) edge_vector_index;
    kv_init(edge_vector_index);
    kvec_t(uint32_t) b;
    kv_init(b);
    uint64_t flag;
    int index;
    
    for (i = 0; i < readNum; i++)  // iterate over all reads
    {
        edge_vector.n = 0;
        edge_vector_index.n = 0;
        for (j = 0; j < (long long)rev_paf[i].length; j++)
        {
            qn = Get_qn(rev_paf[i].buffer[j]);
            tn = Get_tn(rev_paf[i].buffer[j]);
            index = get_specific_overlap(&(paf[tn]), tn, qn);
            if(index != -1)
            {
                flag = tn;
                flag = flag << 32;
                flag = flag | (uint64_t)(index);
                kv_push(uint64_t, edge_vector, flag);
                kv_push(uint64_t, edge_vector_index, j);
            }
        }

        ///based on qn, all edges at edge_vector/edge_vector_index come from different haplotype
        ///but at another direction, all these edges come from the same haplotype
        //here we want to recover these edges
        if((long long)edge_vector_index.n >= rescue_threshold)
        {
            kv_resize(uint32_t, b, edge_vector_index.n);
            b.n = 0;
            for (j = 0; j < (long long)edge_vector_index.n; j++)
            {
                qs = Get_qs(rev_paf[i].buffer[edge_vector_index.a[j]]);
                qe = Get_qe(rev_paf[i].buffer[edge_vector_index.a[j]]);
                kv_push(uint32_t, b, qs<<1);
                kv_push(uint32_t, b, qe<<1|1);
            }

            ks_introsort_uint32_t(b.n, b.a);
            int dp = 0, start = 0, max_dp = 0;
            ma_sub_t max_interval;
            max_interval.s = max_interval.e = 0;
            for (j = 0, dp = 0; j < (long long)b.n; ++j) 
            {
                int old_dp = dp;
                ///if a[j] is qe
                if (b.a[j]&1) 
                {
                    --dp;
                }
                else
                {
                    ++dp;
                } 
                
                /**
                there are two cases: 
                    1. old_dp = dp + 1 (b.a[j] is qe); 2. old_dp = dp - 1 (b.a[j] is qs);
                **/
                ///if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
                if(old_dp < dp) ///b.a[j] is qs
                { 
                    ///case 2, a[j] is qs
                    //here should use dp >= max_dp, instead of dp > max_dp
                    if(dp >= max_dp)
                    {
                        start = b.a[j]>>1;
                        max_dp = dp;
                    }
                } 
                ///else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
                else if (old_dp > dp) ///old_dp > min_dp, b.a[j] is qe
                {
                    if(old_dp == max_dp)
                    {
                        max_interval.s = start;
                        max_interval.e = b.a[j]>>1;
                    }
                }
                // else
                // {
                //     fprintf(stderr, "error\n");
                // }
                
            }

            if(max_dp>= rescue_threshold)
            {
                long long m = 0;
                for (j = 0; j < (long long)edge_vector_index.n; j++)
                {
                    qs = Get_qs(rev_paf[i].buffer[edge_vector_index.a[j]]);
                    qe = Get_qe(rev_paf[i].buffer[edge_vector_index.a[j]]);
                    if(qs <= max_interval.s && qe >= max_interval.e)
                    {
                        edge_vector_index.a[m] = edge_vector_index.a[j];
                        edge_vector.a[m] = edge_vector.a[j];
                        m++;
                    }
                }
                edge_vector_index.n = m;
                edge_vector.n = m;

                ///the read itself do not have these overlaps, but all related reads have
                ///we need to remove all overlaps from rev_paf[i], and then add all overlaps to paf[i]
                remove_overlaps(&(rev_paf[i]), edge_vector_index.a, edge_vector_index.n);
                add_overlaps_from_different_sources(paf, &(paf[i]), edge_vector.a, edge_vector.n);
                revises = revises + edge_vector.n;
            }
        }
    }


    kv_destroy(edge_vector);
    kv_destroy(edge_vector_index);
    kv_destroy(b);
    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] took %0.2fs, rescue edges #: %lld\n\n", 
                                    __func__, Get_T()-startTime, revises);
    }
    
}



long long get_coverage(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, uint64_t n_read)
{
    uint64_t i, j;
    ma_hit_t *h;
    long long R_bases = 0, C_bases = 0;
    for (i = 0; i < n_read; ++i) 
    {
        R_bases += coverage_cut[i].e - coverage_cut[i].s;
        for (j = 0; j < (uint64_t)(sources[i].length); j++)
        {
            h = &(sources[i].buffer[j]);
            C_bases += Get_qe((*h)) - Get_qs((*h));
        }
    }

    return C_bases/R_bases;
}


void pre_clean(ma_hit_t_alloc* sources, ma_sub_t* coverage_cut, asg_t *sg, long long bubble_dist)
{
    while(1)
    {
        int tri_flag = 0;
        ///remove very simple circle
        // (more specifically, loops whose arcs were mostly single edges && also not too long)
        tri_flag += asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, 100, 0);
        ///remove isoloated single read
        // actually this processes small simple bubbles?
        tri_flag += asg_arc_del_single_node_directly(sg, asm_opt.max_short_tip, sources);

        if (!ha_opt_triobin(&asm_opt))  // if not given trio binning files
        {
            tri_flag += asg_arc_del_triangular_advance(sg, bubble_dist);
            ///remove the cross at the bubble carefully, just remove inexact cross
            tri_flag += asg_arc_del_cross_bubble(sg, bubble_dist);
        }
        
        tri_flag += asg_arc_del_single_node_directly(sg, asm_opt.max_short_tip, sources);
        if(tri_flag == 0)
        {
            break;
        }
    }
}


void init_R_to_U(R_to_U* x, uint64_t len)
{
	x->len = len;
	x->index = (uint32_t*)malloc(sizeof(uint32_t)*(x->len));
	memset(x->index, -1, sizeof(uint32_t)*(x->len));
}

void destory_R_to_U(R_to_U* x)
{
	free(x->index);
}

void set_R_to_U(R_to_U* x, uint32_t rID, uint32_t uID, uint32_t is_Unitig)
{
	if(rID >= x->len)
	{
		x->index = (uint32_t*)realloc(x->index, (rID + 1)*sizeof(uint32_t));
        memset(x->index + x->len, -1, sizeof(uint32_t)*((rID + 1) - x->len));
        x->len = rID + 1;
	}

	x->index[rID] = uID & (uint32_t)(0x7fffffff);  // 31bit mask
	x->index[rID] = x->index[rID] | (uint32_t)(is_Unitig<<31);
}


void get_R_to_U(R_to_U* x, uint32_t rID, uint32_t* uID, uint32_t* is_Unitig)
{
	if(rID >= x->len || (x->index[rID] == (uint32_t)(-1)))
	{
		(*uID) = (uint32_t)-1;
		(*is_Unitig) = (uint32_t)-1;
		return;
	}
    
    (*uID) = x->index[rID] & (uint32_t)(0x7fffffff);
	(*is_Unitig) = (x->index[rID]>>31);
}

void transfor_R_to_U(R_to_U* x)
{
    // Update index for all reads, such that 
    // contained reads or alike will have index value
    // pointing to the "longest" parent read.
    // For example, in a nested containment A->B->C 
    // where A is the longest, initial R_to_U might be 
    // B=>A and C=>B. After update, R_to_U should record
    // B=>A and C=>A.
    uint64_t i = 0;
    uint32_t rID, uID, is_Unitig;
    vu32_t san; kv_init(san); kv_resize(uint32_t, san, 16);
    for (i = 0; i < x->len; i++)
    {
        san.n = 0;
        rID = i;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID == (uint32_t)-1) continue;
        if(is_Unitig == 1) continue;


        ///here i/rID is contained in uID
        kv_push(uint32_t, san, rID&0x7fffffff);
        kv_push(uint32_t, san, uID&0x7fffffff);
        rID = uID;
        while (1)
        {
            get_R_to_U(x, rID, &uID, &is_Unitig);
            if(uID == (uint32_t)-1) break;
            if(is_Unitig == 1) break;

            int is_bad = 0;
            for (int j=0; j<san.n; j++){
                if ((uID&0x7fffffff)==san.a[j]){
                    fprintf(stderr, "[E::%s] bad trial of indexing, current uID is %d. dump:\n", 
                            __func__, (int)(uID&0x7fffffff));
                    is_bad = 1;
                }
            }
            if (is_bad){
                for (int j=0; j<san.n; j++)
                    {fprintf(stderr, "rID %d\n", (int)(san.a[j]));}
                break;
            }
            kv_push(uint32_t, san, uID&0x7fffffff);

            rID = uID;
        }

        set_R_to_U(x, i, rID, 0);
    }
    kv_destroy(san);


    /**
    for (i = 0; i < x->len; i++)
    {
        rID = i;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID == (uint32_t)-1) continue;
        if(is_Unitig == 1) continue;
        ///here i/rID is contained in uID
        rID = uID;
        get_R_to_U(x, rID, &uID, &is_Unitig);
        if(uID != (uint32_t)-1)
        {
            fprintf(stderr, "ERROR\n");
        }
    }
    **/
    
}



void asg_delete_node_by_trio(asg_t* sg, uint8_t flag)
{
    uint32_t v, n_vtx = sg->n_seq;
    for (v = 0; v < n_vtx; ++v) 
    {
        if (sg->seq[v].del) continue;
        if (R_INF.trio_flag[v] == AMBIGU) continue;
        if(R_INF.trio_flag[v] != flag) asg_seq_del(sg, v);
    }
    asg_cleanup(sg);
}



void renew_graph_init(ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, asg_t *sg, ma_sub_t *coverage_cut, R_to_U* ruIndex,
uint64_t n_read)
{
    if(sg != NULL) asg_destroy(sg);
    sg = NULL;
    
    if(coverage_cut != NULL) free(coverage_cut);
    coverage_cut = NULL;

    memset(ruIndex->index, -1, sizeof(uint32_t)*(ruIndex->len));


    uint64_t i = 0, j = 0;
    for (i = 0; i < n_read; i++)
    {
        for (j = 0; j < sources[i].length; j++)
        {
            sources[i].buffer[j].del = 0;
        }

        for (j = 0; j < reverse_sources[i].length; j++)
        {
            reverse_sources[i].buffer[j].del = 0;
        }
    }
}




inline uint32_t get_num_edges2existing_nodes_advance(ma_ug_t *ug, asg_t *g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t* oLen,
uint32_t* skip_uId, uint32_t skip_uId_n, uint32_t ignore_trio_flag)
{
    (*oLen) = 0;
    uint32_t qn = query>>1;
    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i, k, occ = 0, uId, is_Unitig, v, w;

    uint32_t trio_flag = R_INF.trio_flag[qn], non_trio_flag = (uint32_t)-1;
    if(ignore_trio_flag == 0)
    {
        if(trio_flag == FATHER) non_trio_flag = MOTHER;
        if(trio_flag == MOTHER) non_trio_flag = FATHER;
    }



    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///sq has already been removed
        ///st must not be removed
        ///g-seq must not be removed
        if(st->del || g->seq[Get_tn(*h)].del) continue;
        if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t.ul>>32) != query) continue;
        get_R_to_U(ruIndex, t.v>>1, &uId, &is_Unitig);
        /****************************may have bugs********************************/
        if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
        /****************************may have bugs********************************/
        for (k = 0; k < skip_uId_n; k++)
        {
            if(uId == skip_uId[k]) break;
        }
        if(k != skip_uId_n) continue;
        // if(uId == skip_uId1) continue;
        // if(uId == skip_uId2) continue;
        /****************************may have bugs********************************/
        if(ug->g->seq[uId].del) continue;
        /****************************may have bugs********************************/
        
        if((t.v != ug->u.a[uId].start) 
            && 
           (t.v != ug->u.a[uId].end))
        {
            continue;
        }
        /**********for debug*************/
        if(t.v == ug->u.a[uId].start)
        {
            ///v = uId<<1;
            v = (uId<<1)^1;
        }

        if(t.v == ug->u.a[uId].end)
        {
            ///v = (uId<<1)^1;
            v = uId<<1;
        }

        if(get_real_length(ug->g, v, NULL)==1)
        {
            get_real_length(ug->g, v, &w);
            if(get_real_length(ug->g, w^1, NULL)==1)
            {
                continue;
            }
        }
        /**********for debug*************/

        occ++;
        (*oLen) += t.ol;
    }

    return occ;
}


inline uint32_t get_edge2existing_node_advance(ma_ug_t *ug, asg_t *g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t* skip_uId, uint32_t skip_uId_n, 
uint32_t* index, asg_arc_t* t, uint32_t ignore_trio_flag)
{
    uint32_t qn = query>>1, uId, is_Unitig, v, w, k;
    int32_t r;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);

    uint32_t trio_flag = R_INF.trio_flag[qn], non_trio_flag = (uint32_t)-1;
    if(ignore_trio_flag == 0)
    {
        if(trio_flag == FATHER) non_trio_flag = MOTHER;
        if(trio_flag == MOTHER) non_trio_flag = FATHER;
    }
    


    for (; (*index) < x->length; (*index)++)
    {
        h = &(x->buffer[(*index)]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///sq has already been removed
        ///st must not be removed
        ///g-seq must not be removed
        if(st->del || g->seq[Get_tn(*h)].del) continue;
        if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t->ul>>32) != query) continue;
        get_R_to_U(ruIndex, t->v>>1, &uId, &is_Unitig);
        /****************************may have bugs********************************/
        if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
        /****************************may have bugs********************************/
        for (k = 0; k < skip_uId_n; k++)
        {
            if(uId == skip_uId[k]) break;
        }
        if(k != skip_uId_n) continue;
        // if(uId == skip_uId1) continue;
        // if(uId == skip_uId2) continue;
        /****************************may have bugs********************************/
        if(ug->g->seq[uId].del) continue;
        /****************************may have bugs********************************/
        
        if((t->v != ug->u.a[uId].start) 
            && 
           (t->v != ug->u.a[uId].end))
        {
            continue;
        }

        /**********for debug*************/
        if(t->v == ug->u.a[uId].start)
        {
            ///v = uId<<1;
            v = (uId<<1)^1;
        }

        if(t->v == ug->u.a[uId].end)
        {
            ///v = (uId<<1)^1;
            v = uId<<1;
        }

        if(get_real_length(ug->g, v, NULL)==1)
        {
            get_real_length(ug->g, v, &w);
            if(get_real_length(ug->g, w^1, NULL)==1)
            {
                continue;
            }
        }
        /**********for debug*************/

        (*index)++;
        return 1;
    }

    return 0;
}



inline uint32_t get_num_edges2existing_nodes(ma_ug_t *ug, asg_t *g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t* oLen,
uint32_t skip_uId1, uint32_t skip_uId2)
{
    (*oLen) = 0;
    uint32_t qn = query>>1;
    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i, occ = 0, uId, is_Unitig, v, w;


    uint32_t trio_flag = R_INF.trio_flag[qn], non_trio_flag = (uint32_t)-1;
    if(trio_flag == FATHER) non_trio_flag = MOTHER;
    if(trio_flag == MOTHER) non_trio_flag = FATHER;


    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///sq has already been removed
        ///st must not be removed
        ///g-seq must not be removed
        if(st->del || g->seq[Get_tn(*h)].del) continue;
        if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t.ul>>32) != query) continue;
        get_R_to_U(ruIndex, t.v>>1, &uId, &is_Unitig);
        /****************************may have bugs********************************/
        if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
        /****************************may have bugs********************************/
        if(uId == skip_uId1) continue;
        if(uId == skip_uId2) continue;
        /****************************may have bugs********************************/
        if(ug->g->seq[uId].del) continue;
        /****************************may have bugs********************************/
        
        if((t.v != ug->u.a[uId].start) 
            && 
           (t.v != ug->u.a[uId].end))
        {
            continue;
        }
        /**********for debug*************/
        if(t.v == ug->u.a[uId].start)
        {
            ///v = uId<<1;
            v = (uId<<1)^1;
        }

        if(t.v == ug->u.a[uId].end)
        {
            ///v = (uId<<1)^1;
            v = uId<<1;
        }

        if(get_real_length(ug->g, v, NULL)==1)
        {
            get_real_length(ug->g, v, &w);
            if(get_real_length(ug->g, w^1, NULL)==1)
            {
                continue;
            }
        }
        /**********for debug*************/

        occ++;
        (*oLen) += t.ol;
    }

    return occ;
}


inline uint32_t get_edge2existing_node(ma_ug_t *ug, asg_t *g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t skip_uId1, uint32_t skip_uId2, 
uint32_t* index, asg_arc_t* t)
{
    uint32_t qn = query>>1, uId, is_Unitig, v, w;
    int32_t r;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);

    uint32_t trio_flag = R_INF.trio_flag[qn], non_trio_flag = (uint32_t)-1;
    if(trio_flag == FATHER) non_trio_flag = MOTHER;
    if(trio_flag == MOTHER) non_trio_flag = FATHER;


    for (; (*index) < x->length; (*index)++)
    {
        h = &(x->buffer[(*index)]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///sq has already been removed
        ///st must not be removed
        ///g-seq must not be removed
        if(st->del || g->seq[Get_tn(*h)].del) continue;
        if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t->ul>>32) != query) continue;
        get_R_to_U(ruIndex, t->v>>1, &uId, &is_Unitig);
        /****************************may have bugs********************************/
        if(uId == (uint32_t)-1 || is_Unitig == 0) continue;
        /****************************may have bugs********************************/
        if(uId == skip_uId1) continue;
        if(uId == skip_uId2) continue;
        /****************************may have bugs********************************/
        if(ug->g->seq[uId].del) continue;
        /****************************may have bugs********************************/
        
        if((t->v != ug->u.a[uId].start) 
            && 
           (t->v != ug->u.a[uId].end))
        {
            continue;
        }

        /**********for debug*************/
        if(t->v == ug->u.a[uId].start)
        {
            ///v = uId<<1;
            v = (uId<<1)^1;
        }

        if(t->v == ug->u.a[uId].end)
        {
            ///v = (uId<<1)^1;
            v = uId<<1;
        }

        if(get_real_length(ug->g, v, NULL)==1)
        {
            get_real_length(ug->g, v, &w);
            if(get_real_length(ug->g, w^1, NULL)==1)
            {
                continue;
            }
        }
        /**********for debug*************/

        (*index)++;
        return 1;
    }

    return 0;
}



uint32_t get_edge_from_source(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query, uint32_t target, asg_arc_t* t)
{
    uint32_t qn = query>>1, i;
    int32_t r;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
    
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t->ul>>32) != query) continue;
        if(t->v != target) continue;
        
        return 1;
    }

    return 0;
}

void append_rId_to_Unitig(asg_t *r_g, ma_ug_t* ug, uint32_t uId, uint32_t rId,
ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    asg_arc_t t;
    t.ul = (uint64_t)-1;
    ma_utg_t *nsu = &(ug->u.a[uId>>1]);
    uint32_t l, i, beg, end;

    if(nsu->m < (nsu->n+1))
    {
        nsu->m = nsu->n+1;
        nsu->a = (uint64_t*)realloc(nsu->a, nsu->m*sizeof(uint64_t));
    }

    if((uId&1) == 0)
    {
        beg = nsu->end^1;
        end = rId;
    }
    
    if((uId&1) == 1)
    {
        beg = rId^1;
        end = nsu->start;
    }

    ///t is the edge that from beg to end
    get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, beg, 
    end, &t);

    ///add rId to the end of nsu
    if((uId&1) == 0)
    {
        nsu->len -= (uint32_t)nsu->a[nsu->n-1];

        l = r_g->seq[rId>>1].len;
        nsu->a[nsu->n] = rId;
        nsu->a[nsu->n] <<=32;
        nsu->a[nsu->n] = nsu->a[nsu->n] | (uint64_t)(l);

        l = asg_arc_len(t);
        nsu->a[nsu->n-1] = nsu->a[nsu->n-1]>>32;
        nsu->a[nsu->n-1] = nsu->a[nsu->n-1]<<32;
        nsu->a[nsu->n-1] = nsu->a[nsu->n-1] | (uint64_t)(l);

        nsu->len = nsu->len + (uint32_t)nsu->a[nsu->n-1] + (uint32_t)nsu->a[nsu->n];
        nsu->end = rId^1;
    
        nsu->n++;
    }

    ///add rId to the start of nsu
    if((uId&1) == 1)
    {
        rId = rId^1;
        for(i = nsu->n; i >= 1; i--)
        {
            nsu->a[i] = nsu->a[i-1];
        }
        
        l = asg_arc_len(t);
        nsu->a[0] = rId;
        nsu->a[0] = nsu->a[0]<<32;
        nsu->a[0] = nsu->a[0] | (uint64_t)(l);
        nsu->len = nsu->len + (uint32_t)nsu->a[0];
        nsu->start = rId;

        nsu->n++;
    }


    set_R_to_U(ruIndex, rId>>1, uId>>1, 1);




    /*************************just for debug**************************/
    uint32_t v, totalLen = 0;
    v = (uint64_t)(nsu->a[0])>>32;
    if(nsu->start != UINT32_MAX && nsu->start != v) fprintf(stderr, "hehe\n");

    v = (uint64_t)(nsu->a[nsu->n-1])>>32;
    if(nsu->end != UINT32_MAX && nsu->end != (v^1)) fprintf(stderr, "haha\n");

    for (i = 0; i < nsu->n; i++)
    {
        l = (uint32_t)(nsu->a[i]);
        totalLen = totalLen + l;
    }
    if(totalLen != nsu->len) fprintf(stderr, "xxxx\n");

    ////fprintf(stderr, "utg%.6dl, dir: %u\n", uId+1, uId&1);
    /*************************just for debug**************************/

}




void lable_all_bubbles(asg_t *r_g, long long bubble_dist)
{
    ///must have this line, otherwise asg_arc_identify_simple_bubbles_multi will be wrong

    asg_cleanup(r_g);
    asg_arc_identify_simple_bubbles_multi(r_g, 0);

    uint32_t v, i, n_vtx = r_g->n_seq * 2;
    buf_t b;
    if (!r_g->is_symm) asg_symm(r_g);
    memset(&b, 0, sizeof(buf_t));

    b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        uint32_t nv = asg_arc_n(r_g, v);
        if(r_g->seq[v>>1].del) continue;
        if(r_g->seq_vis[v] != 0) continue;
        if(nv < 2) continue;
        

        ///if this is a bubble
        ///if(asg_bub_finder_with_del_advance(r_g, v, bubble_dist, &b) == 1)
        if(asg_bub_pop1_primary_trio(r_g, NULL, v, bubble_dist, &b, (uint32_t)-1, (uint32_t)-1, 0))
        {
            //beg is v, end is b.S.a[0]
            //note b.b include end, does not include beg
            for (i = 0; i < b.b.n; i++)
            {
                if(b.b.a[i]==v || b.b.a[i]==b.S.a[0]) continue;

                r_g->seq_vis[b.b.a[i]] = 1;
                r_g->seq_vis[b.b.a[i]^1] = 1;
            }

            r_g->seq_vis[v] = 1;
            r_g->seq_vis[b.S.a[0]^1] = 1;
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
}


void drop_inexact_edegs_at_bubbles(asg_t *r_g, long long bubble_dist)
{
    asg_arc_identify_simple_bubbles_multi(r_g, 0);

    uint32_t v, k, i, n_vtx = r_g->n_seq * 2, nv, flag, n_reduce = 0;
    uint64_t oLen;
    asg_arc_t *av = NULL;
    
    if (!r_g->is_symm) asg_symm(r_g);

    buf_t b;
    memset(&b, 0, sizeof(buf_t));

    kvec_t(uint64_t) e;
    kv_init(e);

    b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
    for (v = 0; v < n_vtx; ++v) 
    {
        nv = asg_arc_n(r_g, v);
        if(r_g->seq[v>>1].del) continue;
        if(r_g->seq_vis[v] != 0) continue;
        if(nv < 2) continue;
        

        ///if this is a bubble
        if(asg_bub_finder_without_del_advance(r_g, v, bubble_dist, &b) == 1)
        {
            //beg is v, end is b.S.a[0]
            //note b.b include end, does not include beg
            for (i = 0; i < b.b.n; i++)
            {
                if(b.b.a[i] == v) continue;
                //note b.b include end, does not include beg
                if(b.b.a[i] == b.S.a[0])
                {
                    b.b.a[i] = b.b.a[i]^1;
                    continue;
                }

                r_g->seq_vis[b.b.a[i]] = 2;
                r_g->seq_vis[b.b.a[i]^1] = 2;
            }

            r_g->seq_vis[v] = 2;
            r_g->seq_vis[b.S.a[0]^1] = 2;
            

            //note b.b include end, does not include beg
            //so push beg into b.b
            kv_push(uint32_t, b.b, v);
            for (i = 0; i < b.b.n; i++)
            {
                nv = asg_arc_n(r_g, b.b.a[i]);
                if(nv <= 1) continue;
                av = asg_arc_a(r_g, b.b.a[i]);
                e.n = 0;
                flag = 0;
                for (k = 0; k < nv; k++)
                {
                    if(av[k].del) continue;
                    if(av[k].el)
                    {
                        flag = 1;
                        continue;
                    }
                    oLen = av[k].ol;
                    oLen = oLen << 32;
                    oLen = oLen | (uint64_t)(k);
                    kv_push(uint64_t, e, oLen);
                }

                if(flag == 0 || e.n == 0) continue;

                if(e.n>1) radix_sort_arch64(e.a, e.a + e.n);

                for (k = 0; k < e.n; k++)
                {
                    if(get_real_length(r_g, av[(uint32_t)e.a[k]].v^1, NULL)<=1) continue;
                    av[(uint32_t)e.a[k]].del = 1;
                    asg_arc_del(r_g, av[(uint32_t)e.a[k]].v^1, av[(uint32_t)e.a[k]].ul>>32^1, 1);
                    n_reduce++;
                }

            }
        }
    }

    free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
    kv_destroy(e);
    if(n_reduce > 0) asg_cleanup(r_g);
}


uint32_t get_corresponding_uId(ma_ug_t *ug, R_to_U* ruIndex, uint32_t t)
{
    uint32_t w, uId, is_Unitig;
    get_R_to_U(ruIndex, t>>1, &uId, &is_Unitig);
    /****************************may have bugs********************************/
    if(uId == ((uint32_t)(-1)) || is_Unitig == 0) return ((uint32_t)(-1));
    if(ug->g->seq[uId].del) return ((uint32_t)(-1));
    /****************************may have bugs********************************/
    w = (uint32_t)-1;
    if(t == ug->u.a[uId].start) w = uId<<1;
    if(t == ug->u.a[uId].end) w = (uId<<1)^1;
    return w;
}

void minor_transitive_reduction(ma_ug_t *ug, R_to_U* ruIndex, asg_arc_t* rbub_edges, uint32_t num)
{
    asg_t* nsg = ug->g;
    uint32_t i, j, k, uId, nv, w;
    asg_arc_t* t = NULL;
    asg_arc_t* p = NULL;
    asg_arc_t *av = NULL;
    ///here all edges from v are saved in rbub_edges
    ///for edges already in graph, need to check del
    ///but for edges in rbub_edges, don't check del
    for (i = 0; i < num; i++)
    {
        t = &rbub_edges[i];
        ///if(t->del) continue;
        uId = get_corresponding_uId(ug, ruIndex, t->v);
        if(uId == ((uint32_t)(-1))) continue;

        nv = asg_arc_n(nsg, uId);
        av = asg_arc_a(nsg, uId);
        for (j = 0; j < nv; j++)
        {
            if(av[j].del) continue;
            w = av[j].v;

            ///note w is the unitig ID
            for (k = 0; k < num; k++)
            {
                p = &rbub_edges[k];
                ///if(p->del) continue;
                ///this line is not necessary at all
                if(k==i) continue;
                ///w is the unitig ID, p->v is the read ID
                if(get_corresponding_uId(ug, ruIndex, p->v) == w) p->del = 1;
            }
        }   
    }    
}




///find contained read with longest overlap
ma_hit_t* get_best_contained_read(ma_ug_t *ug, asg_t *r_g, ma_hit_t_alloc* sources, 
ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query,
uint32_t ignore_trio_flag)
{
    uint32_t qn = query>>1;
    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL, *return_h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(sources[qn]);
    uint32_t i, is_Unitig, contain_rId, contain_uId;
    uint32_t maxOlen = 0;
    asg_t* nsg = ug->g;

    uint32_t trio_flag = R_INF.trio_flag[qn], non_trio_flag = (uint32_t)-1;
    if(ignore_trio_flag == 0)
    {
        if(trio_flag == FATHER) non_trio_flag = MOTHER;
        if(trio_flag == MOTHER) non_trio_flag = FATHER;
    }

    //scan all edges of qn
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///just need deleted edges
        ///sq might be deleted or not
        if(!st->del) continue;
        if(!h->del) continue;
        if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;

        ///tn must be contained in another existing read
        get_R_to_U(ruIndex, Get_tn(*h), &contain_rId, &is_Unitig);
        if(contain_rId == (uint32_t)-1 || is_Unitig == 1) continue;
        if(r_g->seq[contain_rId].del) continue;

        get_R_to_U(ruIndex, contain_rId, &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) continue;
        if(nsg->seq[contain_uId].del) continue;
        ///contain_uId must be a unitig

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t.ul>>32) != query) continue;

        if(t.ol > maxOlen)
        {
            maxOlen = t.ol;
            return_h = h;
        }
    }

    return return_h;
}


int get_contained_reads_chain(ma_hit_t *h, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, 
R_to_U* ruIndex, ma_ug_t *ug, asg_t *r_g, int max_hang, int min_ovlp, uint32_t endRid, uint32_t uId,
kvec_t_u32_warp* chain_buffer, kvec_asg_arc_t_warp* chain_edges, uint32_t* return_ava_cur,
uint32_t* return_ava_ol, uint32_t* return_chainLen, uint32_t thresLen, uint32_t ignore_trio_flag)
{
    (*return_chainLen) = (*return_ava_cur) = (*return_ava_ol) = (uint32_t)-1;
    uint32_t trio_flag, non_trio_flag = (uint32_t)-1, contain_rId, contain_uId, is_Unitig;
    uint32_t chainLen = 0, ava_cur, test_oLen;
    int ql, tl;
    int32_t r;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    asg_arc_t t;
    asg_t* nsg = ug->g;
    chain_buffer->a.n = 0;
    kv_push(uint32_t, chain_buffer->a, uId);
    if(chain_edges) chain_edges->a.n = 0;
    

    ///continue
    ///need to update h, endRid, chainLen, chain_buffer and chain_edges
    while (h)
    {
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        trio_flag = R_INF.trio_flag[Get_qn(*h)];

        ///don't want to edges between different haps
        non_trio_flag = (uint32_t)-1;
        if(ignore_trio_flag == 0)
        {
            if(trio_flag == FATHER) non_trio_flag = MOTHER;
            if(trio_flag == MOTHER) non_trio_flag = FATHER;
            if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) break;
        }

        ///just need deleted edges
        ///sq might be deleted or not
        if(!st->del) break;
        if(!h->del) break;

        ///tn must be contained in another existing read
        get_R_to_U(ruIndex, Get_tn(*h), &contain_rId, &is_Unitig);
        if(contain_rId == (uint32_t)-1 || is_Unitig == 1) break;
        if(r_g->seq[contain_rId].del) break;

        get_R_to_U(ruIndex, contain_rId, &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) break;
        if(nsg->seq[contain_uId].del) break;
        ///contain_uId must be a unitig

        ql = sq->e - sq->s; tl = st->e - st->s;
        r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

        ///if st is contained in sq, or vice verse, skip
        if(r < 0) break;

        ///if sq and v are not in the same direction, skip
        ///endRid is (t.ul>>32), and t.v is a contained read
        if((t.ul>>32) != endRid) break;

        if(chain_edges) kv_push(asg_arc_t, chain_edges->a, t);
        chainLen++;
        kv_push(uint32_t, chain_buffer->a, contain_uId);
        ///endRid is (t.ul>>32), and t.v is a contained read
        ///find edges from t.v to existing unitigs
        ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 
                    ignore_trio_flag);
        //means find an aim
        if(ava_cur > 0)
        {
            /**
            //do nothing if the chainLen == 1
            if(chainLen > 1)
            {
                asg_arc_t t_max;
                ma_hit_t_alloc* x = &(sources[(endRid>>1)]);
                uint32_t ava_ol_max = 0, ava_max = 0, k;
                for (k = 0; k < x->length; k++)
                {
                    h = &(x->buffer[k]);
                    sq = &(coverage_cut[Get_qn(*h)]);
                    st = &(coverage_cut[Get_tn(*h)]);
                    trio_flag = R_INF.trio_flag[Get_qn(*h)];

                    ///don't want to edges between different haps
                    non_trio_flag = (uint32_t)-1;
                    if(trio_flag == FATHER) non_trio_flag = MOTHER;
                    if(trio_flag == MOTHER) non_trio_flag = FATHER;
                    if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;
                    
                    ///just need deleted edges
                    ///sq might be deleted or not
                    if(!st->del) continue;
                    if(!h->del) continue;

                    ///tn must be contained in another existing read
                    get_R_to_U(ruIndex, Get_tn(*h), &contain_rId, &is_Unitig);
                    if(contain_rId == (uint32_t)-1 || is_Unitig == 1) continue;
                    if(r_g->seq[contain_rId].del) continue;

                    get_R_to_U(ruIndex, contain_rId, &contain_uId, &is_Unitig);
                    if(contain_uId == (uint32_t)-1 || is_Unitig != 1) continue;
                    if(nsg->seq[contain_uId].del) continue;
                    ///contain_uId must be a unitig

                    ql = sq->e - sq->s; tl = st->e - st->s;
                    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

                    ///if st is contained in sq, or vice verse, skip
                    if(r < 0) continue;

                    ///if sq and v are not in the same direction, skip
                    if((t.ul>>32) != endRid) continue;

                    ///pop the last contain_uId
                    chain_buffer->a.n--;
                    kv_push(uint32_t, chain_buffer->a, contain_uId);

                    ///note: t.v is a contained read
                    ///we need to find existing reads linked with w
                    ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);

                    if((ava_cur > ava_max) || (ava_cur == ava_max && test_oLen > ava_ol_max))
                    {
                        ava_max = ava_cur;
                        ava_ol_max = test_oLen;
                        t_max = t;
                    }
                }

                if(ava_max > 0)
                {
                    ava_cur = ava_max;
                    test_oLen = ava_ol_max;
                    if(chain_edges)
                    {
                        //pop the last edge
                        chain_edges->a.n--;
                        kv_push(asg_arc_t, chain_edges->a, t_max);
                    } 
                }
                else
                {
                    break;
                }
            }
            **/

            (*return_ava_cur) = ava_cur;
            (*return_ava_ol) = test_oLen;
            (*return_chainLen) = chainLen;
            h = NULL;
            return 1;
        }

        if(chainLen >= thresLen) break;

        ///endRid is (t.ul>>32), and t.v is a contained read
        ///haven't found a existing unitig from t.v
        ///check if t.v can link to a new contained read
        endRid = t.v; 
        h = get_best_contained_read(ug, r_g, sources, coverage_cut, ruIndex, 
            max_hang, min_ovlp, endRid, ignore_trio_flag);
    }

    return 0;    
}

///chainLenThres is used to avoid circle
void rescue_contained_reads_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t chainLenThres, 
uint32_t is_bubble_check, uint32_t is_primary_check, kvec_asg_arc_t_warp* new_rtg_edges, 
kvec_t_u32_warp* new_rtg_nodes)
{
    uint32_t n_vtx, v, k, contain_rId, is_Unitig, uId, rId, endRid, oLen, w, max_oLen, max_oLen_i = (uint32_t)-1;
    uint32_t ava_max, ava_ol_max, ava_min_chain, ava_cur, ava_chainLen, test_oLen, is_update;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    asg_arc_t t, t_max, r_edge;
    t_max.v = t_max.ul = t.v = t.ul = (uint32_t)-1;
    ma_hit_t *h = NULL, *h_max = NULL;
    ma_hit_t_alloc* x = NULL;
    ma_ug_t* ug = NULL;
    uint64_t a_nodes;
    
    kvec_t(asg_arc_t) new_edges;
    kv_init(new_edges);

    kvec_t(asg_arc_t) rbub_edges;
    kv_init(rbub_edges);

    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);

    kvec_t_u32_warp chain_buffer;
    kv_init(chain_buffer.a);

    kvec_asg_arc_t_warp chain_edges;
    kv_init(chain_edges.a);

    if(i_ug != NULL)
    {
        ug = i_ug;
    }
    else
    {
        ug = ma_ug_gen(r_g);
        for (v = 0; v < ug->g->n_seq; v++)
        {
            ug->g->seq[v].c = PRIMARY_LABLE;
        }
    }
    nsg = ug->g;
    
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;
        if(is_primary_check && nsg->seq[v].c==ALTER_LABLE) continue;
        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }




    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        uId = v>>1;
        if(nsg->seq[uId].del) continue;
        if(is_primary_check && nsg->seq[uId].c==ALTER_LABLE) continue;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsu->circ || nsu->start == UINT32_MAX || nsu->end == UINT32_MAX) continue;
        if(get_real_length(nsg, v, NULL) != 0) continue;
        ///we probably don't need this line
        ///if(get_real_length(nsg, v^1, NULL) == 0) continue;
        if(v&1)
        {
            endRid = nsu->start^1;
        }
        else
        {
            endRid = nsu->end^1;
        }

        ///x is the end read of a tip
        ///find all overlap of x
        x = &(sources[(endRid>>1)]);
        ava_ol_max = ava_max = 0; ava_min_chain = (uint32_t)-1;
        h_max = NULL;
        for (k = 0; k < x->length; k++)
        {
            ///h is the edge of endRid
            h = &(x->buffer[k]);
            ///means we found a contained read
            if(get_contained_reads_chain(h, sources, coverage_cut, ruIndex, ug, r_g, 
            max_hang, min_ovlp, endRid, uId, &chain_buffer, NULL, &ava_cur, &test_oLen, 
            &ava_chainLen, chainLenThres, 1))
            {
                is_update = 0;
                if(ava_chainLen < ava_min_chain)
                {
                    is_update = 1;
                }
                else if(ava_chainLen == ava_min_chain)
                {
                    if(ava_cur > ava_max)
                    {
                        is_update = 1;
                    }
                    else if(ava_cur == ava_max && test_oLen > ava_ol_max)
                    {
                        is_update = 1;
                    }
                }
                if(is_update)
                {
                    ava_min_chain = ava_chainLen;
                    ava_max = ava_cur;
                    ava_ol_max = test_oLen;
                    h_max = h;
                }
            }
        }





        if(ava_max > 0)
        {
            //get all edges between contained reads
            get_contained_reads_chain(h_max, sources, coverage_cut, ruIndex, ug, r_g, 
            max_hang, min_ovlp, endRid, uId, &chain_buffer, &chain_edges, &ava_cur, &test_oLen, 
            &ava_chainLen, chainLenThres, 1);
            if(chain_edges.a.n < 1) continue;
            ///the last cantained read
            t_max = chain_edges.a.a[chain_edges.a.n-1];

            k = 0; rbub_edges.n = 0;
            ///edges from the last contained read to other unitigs
            while(get_edge2existing_node_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t_max.v, chain_buffer.a.a, chain_buffer.a.n, &k, &r_edge, 1))
            {
               kv_push(asg_arc_t, rbub_edges, r_edge);
            }

            ///need to do transitive reduction
            ///note here is different to standard transitive reduction
            minor_transitive_reduction(ug, ruIndex, rbub_edges.a, rbub_edges.n);
            if(is_primary_check)
            {
                max_oLen = 0; max_oLen_i = (uint32_t)-1;
                for (k = 0; k < rbub_edges.n; k++)
                {
                    t = rbub_edges.a[k];
                    if(t.del) continue;
                    if(t.ol > max_oLen)
                    {
                        max_oLen = t.ol;
                        max_oLen_i = k;
                    }
                }

                if(max_oLen_i == (uint32_t)-1) continue;
                t = rbub_edges.a[max_oLen_i];
                w = get_corresponding_uId(ug, ruIndex, t.v);
                if(w == ((uint32_t)(-1))) continue;
                if(get_real_length(nsg, w^1, NULL)!=0) continue;
            }
            
            //recover all contained reads
            for (k = 0; k < chain_edges.a.n; k++)
            {
                t_max = chain_edges.a.a[k];
                ///save all infor for reverting
                get_R_to_U(ruIndex, t_max.v>>1, &contain_rId, &is_Unitig);
                a_nodes=contain_rId; 
                a_nodes=a_nodes<<32; 
                a_nodes=a_nodes|((uint64_t)(t_max.v>>1));
                kv_push(uint64_t, u_vecs.a, a_nodes);


                //(t_max.ul>>32)----->(t_max.v/r_edge.ul>>32)----->r_edge.v
                //need to recover (t_max.v/r_edge.ul>>32), 1. set .del = 0 2. set ruIndex
                r_g->seq[t_max.v>>1].del = 0;
                coverage_cut[t_max.v>>1].del = 0;
                coverage_cut[t_max.v>>1].c = PRIMARY_LABLE;
                //add recover to the end of the unitig
                append_rId_to_Unitig(r_g, ug, v, t_max.v, sources, coverage_cut, ruIndex, max_hang, 
                min_ovlp);

                get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t_max.ul>>32), t_max.v, &t);
                kv_push(asg_arc_t, new_edges, t);

                get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t_max.v^1), ((t_max.ul>>32)^1), &t);
                kv_push(asg_arc_t, new_edges, t);
            }

            if(is_primary_check)
            {
                t = rbub_edges.a[max_oLen_i];
                w = get_corresponding_uId(ug, ruIndex, t.v);
                if(w == ((uint32_t)(-1))) continue;
                kv_push(asg_arc_t, new_edges, t);
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);

                get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t.v^1), ((t.ul>>32)^1), &t);
                kv_push(asg_arc_t, new_edges, t);
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
            }
            else
            {
                for (k = 0; k < rbub_edges.n; k++)
                {
                    t = rbub_edges.a[k];
                    if(t.del) continue;
                    w = get_corresponding_uId(ug, ruIndex, t.v);
                    if(w == ((uint32_t)(-1))) continue;
            
                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);


                    get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                    (t.v^1), ((t.ul>>32)^1), &t);
                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
                }
            }
            
        }
    }

    asg_arc_t* p = NULL;
    if(is_bubble_check)
    {
        
        for (k = 0; k < new_edges.n; k++)
        {
            p = asg_arc_pushp(r_g);
            *p = new_edges.a[k];
        }

        if(new_edges.n != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
            asg_cleanup(r_g);
        }
        
        while(1)
        {
            int tri_flag = 0;
            ///remove very simple circle
            tri_flag += asg_arc_del_simple_circle_untig(sources, coverage_cut, r_g, 100, 0);
            if (!ha_opt_triobin(&asm_opt))
            {
                ///remove isoloated single read
                tri_flag += asg_arc_del_triangular_advance(r_g, bubble_dist);
                ///remove the cross at the bubble carefully, just remove inexact cross
                tri_flag += asg_arc_del_cross_bubble(r_g, bubble_dist);
            }
            if(tri_flag == 0)
            {
                break;
            }
        }

        lable_all_bubbles(r_g, bubble_dist);


        for (k = 0; k < new_edges.n; k++)
        {
            v = new_edges.a[k].ul>>32;
            w = new_edges.a[k].v;

            if(r_g->seq[v>>1].del) continue;
            if(r_g->seq[w>>1].del) continue;
            ///if this edge is at a bubble
            if(r_g->seq_vis[v]!=0 && r_g->seq_vis[w^1]!=0) continue;

            asg_arc_del(r_g, v, w, 1);
            asg_arc_del(r_g, w^1, v^1, 1);
        }
        asg_cleanup(r_g);


        // w=max_contain_rId; w=w<<32; w=w|(t_max.v>>1);
        // kv_push(uint64_t, u_vecs.a, w);
        for (k = 0; k < u_vecs.a.n; k++)
        {
            w = (uint32_t)u_vecs.a.a[k];
            w = w<<1;
            if(!r_g->seq[w>>1].del && get_real_length(r_g, w, NULL)!=0) continue;
            if(!r_g->seq[w>>1].del && get_real_length(r_g, (w^1), NULL)!=0) continue;
            w=w>>1;
            r_g->seq[w].del = 1;
            coverage_cut[w].del = 1;
            set_R_to_U(ruIndex, ((uint32_t)(u_vecs.a.a[k])), (u_vecs.a.a[k]>>32), 0);
        }
    }

    ///actually we don't update sources
    ///during the the primary_check, we need to recover everything for r_g, coverage_cut and ruIndex
    if(is_primary_check)
    {
        ///don't remove nodes and edges
        /**
        for (k = 0; k < u_vecs.a.n; k++)
        {
            w = (uint32_t)u_vecs.a.a[k];
            r_g->seq[w].del = 1;
            coverage_cut[w].del = 1;
            set_R_to_U(ruIndex, ((uint32_t)(u_vecs.a.a[k])), (u_vecs.a.a[k]>>32), 0);
        }
        **/
       if(new_rtg_nodes)
       {
           for(k = 0; k < u_vecs.a.n; k++)
           {
                w = (uint32_t)u_vecs.a.a[k];
                kv_push(uint32_t, new_rtg_nodes->a, w);
           }
       }

       
        for (k = 0; k < new_edges.n; k++)
        {
            p = asg_arc_pushp(r_g);
            *p = new_edges.a[k];
            if(new_rtg_edges) kv_push(asg_arc_t, new_rtg_edges->a, new_edges.a[k]);
        }

        if(new_edges.n != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
            asg_cleanup(r_g);
        }
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }


    if(i_ug == NULL) ma_ug_destroy(ug);
    kv_destroy(new_edges);
    kv_destroy(rbub_edges);
    kv_destroy(u_vecs.a);
    kv_destroy(chain_buffer.a);
    kv_destroy(chain_edges.a);
    
}



void rescue_missing_overlaps_aggressive(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t is_bubble_check, uint32_t is_primary_check,
kvec_asg_arc_t_warp* new_rtg_edges)
{
    uint32_t n_vtx, v, k, is_Unitig, uId, rId, endRid, oLen, w, max_oLen, max_oLen_i;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    ma_ug_t *ug = NULL;
    asg_arc_t t, r_edge;
    kvec_t(asg_arc_t) new_edges;
    kv_init(new_edges);

    kvec_t(asg_arc_t) rbub_edges;
    kv_init(rbub_edges);

    if(i_ug != NULL)
    {
        ug = i_ug;
    }
    else
    {
        ug = ma_ug_gen(r_g);
        for (v = 0; v < ug->g->n_seq; v++)
        {
            ug->g->seq[v].c = PRIMARY_LABLE;
        }
    }
    nsg = ug->g;
    
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;
        if(is_primary_check && nsg->seq[v].c==ALTER_LABLE) continue;

        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }


    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        uId = v>>1;
        if(nsg->seq[uId].del) continue;
        if(is_primary_check && nsg->seq[uId].c == ALTER_LABLE) continue;

        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsu->circ || nsu->start == UINT32_MAX || nsu->end == UINT32_MAX) continue;
        if(get_real_length(nsg, v, NULL) != 0) continue;
        
        if(v&1)
        {
            endRid = nsu->start^1;
        }
        else
        {
            endRid = nsu->end^1;
        }


        k = 0;
        rbub_edges.n = 0;
        while(get_edge2existing_node_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                max_hang, min_ovlp, endRid, &uId, 1, &k, &r_edge, 1))
        {
            kv_push(asg_arc_t, rbub_edges, r_edge);
        }



        if(rbub_edges.n > 0)
        {
            ///need to do transitive reduction
            ///note here is different to standard transitive reduction
            minor_transitive_reduction(ug, ruIndex, rbub_edges.a, rbub_edges.n);
            if(is_primary_check)
            {
                max_oLen = 0; max_oLen_i = (uint32_t)-1;
                for (k = 0; k < rbub_edges.n; k++)
                {
                    t = rbub_edges.a[k];
                    if(t.del) continue;
                    if(t.ol > max_oLen)
                    {
                        max_oLen = t.ol;
                        max_oLen_i = k;
                    }
                }

                if(max_oLen_i == (uint32_t)-1) continue;
                t = rbub_edges.a[max_oLen_i];
                w = get_corresponding_uId(ug, ruIndex, t.v);
                if(w == ((uint32_t)(-1))) continue;
                if(get_real_length(nsg, w^1, NULL)!=0) continue;

                kv_push(asg_arc_t, new_edges, t);
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);


                get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t.v^1), ((t.ul>>32)^1), &t);
                kv_push(asg_arc_t, new_edges, t);
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
            }
            else
            {
                for (k = 0; k < rbub_edges.n; k++)
                {
                    t = rbub_edges.a[k];
                    if(t.del) continue;
                    w = get_corresponding_uId(ug, ruIndex, t.v);
                    if(w == ((uint32_t)(-1))) continue;
                    // get_R_to_U(ruIndex, t.v>>1, &uId, &is_Unitig);
                    // w = (uint32_t)-1;
                    // if(t.v == ug->u.a[uId].start) w = uId<<1;
                    // if(t.v == ug->u.a[uId].end) w = (uId<<1)^1;
                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);


                    get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                    (t.v^1), ((t.ul>>32)^1), &t);
                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
                }
            }

        }
    }

    asg_arc_t* p = NULL;
    if(is_bubble_check)
    {
        
        for (k = 0; k < new_edges.n; k++)
        {
            p = asg_arc_pushp(r_g);
            *p = new_edges.a[k];
        }

        if(new_edges.n != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
            asg_cleanup(r_g);
        }
        

        while(1)
        {
            int tri_flag = 0;
            ///remove very simple circle
            tri_flag += asg_arc_del_simple_circle_untig(sources, coverage_cut, r_g, 100, 0);
            if (!ha_opt_triobin(&asm_opt))
            {
                ///remove isoloated single read
                tri_flag += asg_arc_del_triangular_advance(r_g, bubble_dist);
                ///remove the cross at the bubble carefully, just remove inexact cross
                tri_flag += asg_arc_del_cross_bubble(r_g, bubble_dist);
            }
            if(tri_flag == 0)
            {
                break;
            }
        }

        lable_all_bubbles(r_g, bubble_dist);


        for (k = 0; k < new_edges.n; k++)
        {
            v = new_edges.a[k].ul>>32;
            w = new_edges.a[k].v;

            if(r_g->seq[v>>1].del) continue;
            if(r_g->seq[w>>1].del) continue;
            ///if this edge is at a bubble
            if(r_g->seq_vis[v]!=0 && r_g->seq_vis[w^1]!=0) continue;

            asg_arc_del(r_g, v, w, 1);
            asg_arc_del(r_g, w^1, v^1, 1);
        }
        asg_cleanup(r_g);
    }

    if(is_primary_check)
    {
        for (k = 0; k < new_edges.n; k++)
        {
            p = asg_arc_pushp(r_g);
            *p = new_edges.a[k];
            if(new_rtg_edges) kv_push(asg_arc_t, new_rtg_edges->a, new_edges.a[k]);
        }

        if(new_edges.n != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
            asg_cleanup(r_g);
        }
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    if(i_ug == NULL) ma_ug_destroy(ug);

    kv_destroy(new_edges);
    kv_destroy(rbub_edges);
    /*************************just for debug**************************/
    // debug_utg_graph(ug, r_g, 0, 0);
    /*************************just for debug**************************/
}

void update_unitig(long long step, long long init, ma_utg_t* nsu, asg_t *r_g, 
kvec_asg_arc_t_warp* recover_edges, uint32_t update_mode)
{
    uint64_t l, k;
    uint32_t v;
    ///recover
    if(update_mode == 0)
    {
        if(step == 1) nsu->start = ((uint64_t)(nsu->a[0])>>32);
        if(step == -1) nsu->end = ((uint64_t)(nsu->a[nsu->n-1])>>32)^1;
        if(recover_edges)
        {
            ///the first node has not been deleted, others has been deleted
            for (init = init - step; init >= 0 && init < (long long)nsu->n; init = init - step)
            {
                v = (uint64_t)nsu->a[init]>>32;
                r_g->seq[v>>1].del = 0;
                for (k = 0; k < recover_edges->a.n; k++)
                {
                    if(((v>>1)==(recover_edges->a.a[k].ul>>33)) || 
                                        ((v>>1)==(recover_edges->a.a[k].v>>1)))
                    {
                        recover_edges->a.a[k].del = 0;
                    }
                }
            }
        }
    }
    else ///update
    {
        //end of unitig
        if(step == -1 && init != (long long)((long long)nsu->n - 1))
        {
            l = r_g->seq[(nsu->a[init]>>33)].len;
            nsu->n = init+1;
            nsu->a[nsu->n-1] = nsu->a[nsu->n-1]>>32;
            nsu->a[nsu->n-1] = nsu->a[nsu->n-1]<<32;
            nsu->a[nsu->n-1] = nsu->a[nsu->n-1] | (uint64_t)(l);
        }

        //beg of unitig
        if(step == 1 && init != 0)
        {
            for(k = init; k < nsu->n; k++)
            {
                nsu->a[k-init] = nsu->a[k];
            }
            nsu->n = nsu->n - init;
        }

        nsu->len = 0;
        for(k = 0; k < nsu->n; k++)
        {
            nsu->len += (uint32_t)nsu->a[k];
        }
    }
}
void rescue_missing_overlaps_backward(ma_ug_t *i_ug, asg_t *r_g, ma_hit_t_alloc* sources, ma_sub_t *coverage_cut,
R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist, uint32_t backward_steps, 
uint32_t is_bubble_check, uint32_t is_primary_check)
{
    uint32_t v, vId, dir, k, cur_backward_steps, round, is_Unitig, uId, rId, endRid, oLen, w, max_oLen, max_oLen_i, mode, nv;
    uint64_t tmp;
    long long init, step = 0;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    ma_ug_t *ug = NULL;
    asg_arc_t *av = NULL;
    asg_arc_t t, r_edge;
    kvec_t(asg_arc_t) new_edges;
    kv_init(new_edges);

    kvec_asg_arc_t_warp recover_edges;
    kv_init(recover_edges.a);

    kvec_t(asg_arc_t) rbub_edges;
    kv_init(rbub_edges);

    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);

    if(i_ug != NULL)
    {
        ug = i_ug;
    }
    else
    {
        ug = ma_ug_gen(r_g);
        for (v = 0; v < ug->g->n_seq; v++)
        {
            ug->g->seq[v].c = PRIMARY_LABLE;
        }
    }
    nsg = ug->g;
    
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;
        if(is_primary_check && nsg->seq[v].c==ALTER_LABLE) continue;

        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }



    for (vId = 0; vId < nsg->n_seq; vId++)
    {
        rbub_edges.n = round = 0;
        for (dir = 0; dir < 2; dir++)
        {
            if(rbub_edges.n > 0)
            {
                cur_backward_steps = nsu->n - round - 1;
                if(cur_backward_steps > backward_steps)
                {
                    cur_backward_steps = backward_steps;
                }
            }
            else
            {
                cur_backward_steps = backward_steps;
            }
            
            

            v = vId; v = v<<1; v = v | dir;
            uId = v>>1;
            if(nsg->seq[uId].del) continue;
            if(is_primary_check && nsg->seq[uId].c == ALTER_LABLE) continue;

            nsu = &(ug->u.a[uId]);
            if(nsu->m == 0) continue;
            if(nsu->circ || nsu->start == UINT32_MAX || nsu->end == UINT32_MAX) continue;
            if(get_real_length(nsg, v, NULL) != 0) continue;
            ////if(get_real_length(nsg, v^1, NULL) == 0) continue;
            ///that means this unitig has been changed
            if(nsu->start!=((uint64_t)(nsu->a[0])>>32)) continue;
            if((nsu->end^1)!=((uint64_t)(nsu->a[nsu->n-1])>>32)) continue;
            
            if(v&1)
            {
                init = 0;
                step = 1;
                mode = 1;
                //endRid = nsu->start^1;
            }
            else
            {
                init = nsu->n - 1;
                step = -1;
                mode = 0;
                ///endRid = nsu->end^1;
            }


            rbub_edges.n = 0;
            for (round = 0; round < cur_backward_steps && init >= 0 && init < (long long)nsu->n; 
                                                                    init = init + step, round++)
            {
                endRid = ((uint64_t)(nsu->a[init]))>>32;
                endRid = endRid^mode;

                k = 0;
                rbub_edges.n = 0;
                while(get_edge2existing_node_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, endRid, &uId, 1, &k, &r_edge, 1))
                {
                    kv_push(asg_arc_t, rbub_edges, r_edge);
                }
                if(rbub_edges.n > 0) break;
            }

            if(rbub_edges.n > 0)
            {
                if(mode == 1) nsu->start = endRid^1;
                if(mode == 0) nsu->end = endRid^1;
                //save for revert
                tmp = mode; tmp = tmp <<31; tmp = tmp | (uint64_t)(init); tmp = tmp << 32; tmp = tmp | uId;
                kv_push(uint64_t, u_vecs.a, tmp);



                ///need to do transitive reduction
                ///note here is different to standard transitive reduction
                minor_transitive_reduction(ug, ruIndex, rbub_edges.a, rbub_edges.n);
                if(is_primary_check)
                {
                    max_oLen = 0; max_oLen_i = (uint32_t)-1;
                    for (k = 0; k < rbub_edges.n; k++)
                    {
                        t = rbub_edges.a[k];
                        if(t.del) continue;
                        if(t.ol > max_oLen)
                        {
                            max_oLen = t.ol;
                            max_oLen_i = k;
                        }
                    }

                    if(max_oLen_i == (uint32_t)-1) continue;
                    t = rbub_edges.a[max_oLen_i];
                    w = get_corresponding_uId(ug, ruIndex, t.v);
                    if(w == ((uint32_t)(-1))) continue;
                    if(get_real_length(nsg, w^1, NULL)!=0) continue;

                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);


                    get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                    (t.v^1), ((t.ul>>32)^1), &t);
                    kv_push(asg_arc_t, new_edges, t);
                    oLen = t.ol;
                    asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
                }
                else
                {
                    ///modify read graph
                    for (init = init - step; init >= 0 && init < (long long)nsu->n; init = init - step)
                    {
                        w = ((uint64_t)(nsu->a[init]))>>32;
                        nv = asg_arc_n(r_g, w);
                        av = asg_arc_a(r_g, w);
                        for (k = 0; k < nv; k++)
                        {
                            if(av[k].del) continue;
                            kv_push(asg_arc_t, recover_edges.a, av[k]);
                            if(asg_get_arc(r_g, av[k].v^1, av[k].ul>>32^1, &t)==0)
                            {
                                fprintf(stderr, "[E::%s]error\n", __func__);
                            } 
                            kv_push(asg_arc_t, recover_edges.a, t);
                        }
                        

                        nv = asg_arc_n(r_g, w^1);
                        av = asg_arc_a(r_g, w^1);
                        for (k = 0; k < nv; k++)
                        {
                            if(av[k].del) continue;
                            kv_push(asg_arc_t, recover_edges.a, av[k]);
                            if(asg_get_arc(r_g, av[k].v^1, av[k].ul>>32^1, &t)==0)
                            {
                                fprintf(stderr, "[E::%s]error\n", __func__);
                            } 
                            kv_push(asg_arc_t, recover_edges.a, t);
                        }

                        ///w = ((uint64_t)(nsu->a[init]))>>32;
                        asg_seq_del(r_g, w>>1);
                    }


                    
                    for (k = 0; k < rbub_edges.n; k++)
                    {
                        t = rbub_edges.a[k];
                        if(t.del) continue;
                        w = get_corresponding_uId(ug, ruIndex, t.v);
                        if(w == ((uint32_t)(-1))) continue;
                        // get_R_to_U(ruIndex, t.v>>1, &uId, &is_Unitig);
                        // w = (uint32_t)-1;
                        // if(t.v == ug->u.a[uId].start) w = uId<<1;
                        // if(t.v == ug->u.a[uId].end) w = (uId<<1)^1;
                        kv_push(asg_arc_t, new_edges, t);
                        oLen = t.ol;
                        asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);


                        get_edge_from_source(sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                        (t.v^1), ((t.ul>>32)^1), &t);
                        kv_push(asg_arc_t, new_edges, t);
                        oLen = t.ol;
                        asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);
                    }
                }

            }
        }
    }
    
    
    if(is_bubble_check)
    {
        asg_arc_t* p = NULL;
        for (k = 0; k < new_edges.n; k++)
        {
            p = asg_arc_pushp(r_g);
            *p = new_edges.a[k];
        }

        if(new_edges.n != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
            asg_cleanup(r_g);
        }
        
        while(1)
        {
            int tri_flag = 0;
            ///remove very simple circle
            tri_flag += asg_arc_del_simple_circle_untig(sources, coverage_cut, r_g, 100, 0);
            if (!ha_opt_triobin(&asm_opt))
            {
                ///remove isoloated single read
                tri_flag += asg_arc_del_triangular_advance(r_g, bubble_dist);
                ///remove the cross at the bubble carefully, just remove inexact cross
                tri_flag += asg_arc_del_cross_bubble(r_g, bubble_dist);
            }
            if(tri_flag == 0)
            {
                break;
            }
        }

        lable_all_bubbles(r_g, bubble_dist);


        for (k = 0; k < new_edges.n; k++)
        {
            v = new_edges.a[k].ul>>32;
            w = new_edges.a[k].v;

            if(r_g->seq[v>>1].del) continue;
            if(r_g->seq[w>>1].del) continue;
            ///if this edge is at a bubble
            if(r_g->seq_vis[v]!=0 && r_g->seq_vis[w^1]!=0) continue;

            asg_arc_del(r_g, v, w, 1);
            asg_arc_del(r_g, w^1, v^1, 1);
        }

        asg_cleanup(r_g);
        

        for (k = 0; k < recover_edges.a.n; k++)
        {
            recover_edges.a.a[k].del = 1;
        }

        for (k = 0; k < u_vecs.a.n; k++)
        {
            mode = (uint64_t)u_vecs.a.a[k]>>63;
            if(mode == 1) step = 1;
            if(mode == 0) step = -1;
            init = (uint64_t)((uint64_t)u_vecs.a.a[k]>>32)&((uint64_t)(0x7fffffff));
            uId = (uint32_t)u_vecs.a.a[k];
            nsu = &(ug->u.a[uId]);

            endRid = ((uint64_t)(nsu->a[init]))>>32;
            endRid = endRid^mode;
            /********************for debug**********************/
            // fprintf(stderr, "n: %u, k: %u, mode: %u, init: %lld, step: %lld, uId: %u, endRid: %u\n", 
            // (uint32_t)u_vecs.a.n, k, mode, init, step, uId, endRid);
            // fprintf(stderr, "nsu->start: %u, nsu->end: %u\n", 
            // nsu->start, nsu->end);
            // if(mode == 1 && nsu->start != (endRid^1)) fprintf(stderr, "ERROR\n");
            // if(mode == 0 && nsu->end != (endRid^1)) fprintf(stderr, "ERROR\n");
            // update_unitig(step, init, nsu, r_g, &recover_edges, 0);
            /********************for debug**********************/
            if(get_real_length(r_g, endRid, NULL) > 0)
            {
                update_unitig(step, init, nsu, r_g, &recover_edges, 1);
                // fprintf(stderr, "endRid: %u, %.*s, uId: %u\n", 
                // endRid>>1, (int)Get_NAME_LENGTH((R_INF), v>>1), Get_NAME((R_INF), v>>1), uId);
            }
            else
            {
                update_unitig(step, init, nsu, r_g, &recover_edges, 0);
            }
        }

        uint64_t recov_occ = 0;
        for (k = 0; k < recover_edges.a.n; k++)
        {
            if(recover_edges.a.a[k].del) continue;
            p = asg_arc_pushp(r_g);
            *p = recover_edges.a.a[k];
            recov_occ++;
        }

        if(recov_occ != 0)
        {
            free(r_g->idx);
            r_g->idx = 0;
            r_g->is_srt = 0;
        }

        asg_cleanup(r_g);
    }

    if(is_primary_check)
    {
        for (k = 0; k < u_vecs.a.n; k++)
        {
            mode = (uint64_t)u_vecs.a.a[k]>>63;
            if(mode == 1) step = 1;
            if(mode == 0) step = -1;
            init = (uint64_t)((uint64_t)u_vecs.a.a[k]>>32)&((uint64_t)(0x7fffffff));
            uId = (uint32_t)u_vecs.a.a[k];
            nsu = &(ug->u.a[uId]);

            endRid = ((uint64_t)(nsu->a[init]))>>32;
            endRid = endRid^mode;

            update_unitig(step, init, nsu, r_g, &recover_edges, 1);
        }
    }


    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    if(i_ug == NULL) ma_ug_destroy(ug);

    kv_destroy(new_edges);
    kv_destroy(recover_edges.a);
    kv_destroy(rbub_edges);
    kv_destroy(u_vecs.a);
}




void rescue_wrong_overlaps_to_unitigs(ma_ug_t *i_ug, asg_t *r_g,  ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, long long bubble_dist,
kvec_asg_arc_t_warp* keep_edges)
{
    uint32_t n_vtx, v, k, is_Unitig, uId, rId, endRid, oLen, w;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    ma_ug_t *ug = NULL;
    asg_arc_t t, r_edge;
    ma_hit_t_alloc* x = NULL;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    int32_t r;

    kvec_t(asg_arc_t) new_utg_edges;
    kv_init(new_utg_edges);

    kvec_t(asg_arc_t) new_rtg_edges;
    kv_init(new_rtg_edges);

    kvec_t(asg_arc_t) rbub_edges;
    kv_init(rbub_edges);

    if(i_ug != NULL)
    {
        ug = i_ug;
    }
    else
    {
        ug = ma_ug_gen(r_g);
        for (v = 0; v < ug->g->n_seq; v++)
        {
            ug->g->seq[v].c = PRIMARY_LABLE;
        }
    }

    nsg = ug->g;

    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;

        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }

    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        uId = v>>1;
        if(nsg->seq[uId].del) continue;

        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsu->circ || nsu->start == UINT32_MAX || nsu->end == UINT32_MAX) continue;
        if(get_real_length(nsg, v, NULL) != 0) continue;
        


        if(v&1)
        {
            endRid = nsu->start^1;
        }
        else
        {
            endRid = nsu->end^1;
        }



        /****************************may have bugs********************************/
        x = &(sources[(endRid>>1)]);
        for (k = 0; k < x->length; k++)
        {
            ///h is the edge of endRid
            h = &(x->buffer[k]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            ///sq, st, and h cannot be deleted
            if(sq->del || st->del || h->del) continue;
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                            asm_opt.max_hang_rate, min_ovlp, &t);
            if(r < 0) continue;
            if((t.ul>>32) != endRid) continue;
            break;
        }
        ///means there is >0 edges
        if(k != x->length) continue;
        /****************************may have bugs********************************/




        k = 0;
        rbub_edges.n = 0;
        ///note here use reverse_sources instead of sources
        while(get_edge2existing_node_advance(ug, r_g, reverse_sources, coverage_cut, ruIndex, 
                max_hang, min_ovlp, endRid, &uId, 1, &k, &r_edge, 1))
        {
            kv_push(asg_arc_t, rbub_edges, r_edge);
        }


        
        if(rbub_edges.n > 0)
        {
            ///need to do transitive reduction
            ///note here is different to standard transitive reduction
            minor_transitive_reduction(ug, ruIndex, rbub_edges.a, rbub_edges.n);

            for (k = 0; k < rbub_edges.n; k++)
            {
                t = rbub_edges.a[k];
                if(t.del) continue;
                w = get_corresponding_uId(ug, ruIndex, t.v);
                if(w == ((uint32_t)(-1))) continue;
                ///if(uId == 2474) fprintf(stderr, "###uId: %u, v: %u, r_edge.v>>1: %u\n", uId, v, r_edge.v>>1);
                kv_push(asg_arc_t, new_rtg_edges, t);
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, oLen, t.strong, t.el, t.no_l_indel);
                

                get_edge_from_source(reverse_sources, coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t.v^1), ((t.ul>>32)^1), &t);
                kv_push(asg_arc_t, new_rtg_edges, t);
                ///if edge does not exist in reverse, oLen won't change, it is what we want
                oLen = t.ol;
                asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, oLen, t.strong, t.el, t.no_l_indel);

                t.ul = v; t.ul = t.ul<<32; t.v = w;
                kv_push(asg_arc_t, new_utg_edges, t);


                t.ul = w^1; t.ul = t.ul<<32; t.v = v^1;
                kv_push(asg_arc_t, new_utg_edges, t);
            }
        }
        
    }


    nsg->seq_vis = (uint8_t*)calloc(nsg->n_seq*2, sizeof(uint8_t));    
    lable_all_bubbles(nsg, bubble_dist);

    /*********************************for debug**************************************/
    // uint32_t v_uId, w_uId;
    // if(new_utg_edges.n != new_rtg_edges.n) fprintf(stderr, "ERROR 1\n");
    // for (k = 0; k < new_utg_edges.n; k++)
    // {
    //     v = new_rtg_edges.a[k].ul>>32;
    //     v_uId = (get_corresponding_uId(ug, ruIndex, v^1)^1);
    //     w = new_rtg_edges.a[k].v;
    //     w_uId = get_corresponding_uId(ug, ruIndex, w);

    //     v = new_utg_edges.a[k].ul>>32;
    //     w = new_utg_edges.a[k].v;

    //     if(v!=v_uId || w!= w_uId)
    //     {
    //         fprintf(stderr, "\n(%u) ERROR 2\n", k);
    //         fprintf(stderr, "v>>1: %u, v&1: %u, ug->u.a[v>>1].n: %u\n", 
    //                                                     v>>1, v&1, ug->u.a[v>>1].n);
    //         fprintf(stderr, "ug->u.a[v>>1].start>>1: %u, ug->u.a[v>>1].start&1: %u\n", 
    //                                         ug->u.a[v>>1].start>>1, ug->u.a[v>>1].start&1);
    //         fprintf(stderr, "ug->u.a[v>>1].end>>1: %u, ug->u.a[v>>1].end&1: %u\n", 
    //                                         ug->u.a[v>>1].end>>1, ug->u.a[v>>1].end&1);
    //         fprintf(stderr, "w>>1: %u, w&1: %u, ug->u.a[w>>1].n: %u\n", 
    //                                                     w>>1, w&1, ug->u.a[w>>1].n);
    //         fprintf(stderr, "ug->u.a[w>>1].start>>1: %u, ug->u.a[w>>1].start&1: %u\n", 
    //                                         ug->u.a[w>>1].start>>1, ug->u.a[w>>1].start&1);
    //         fprintf(stderr, "ug->u.a[w>>1].end>>1: %u, ug->u.a[w>>1].end&1: %u\n", 
    //                                         ug->u.a[w>>1].end>>1, ug->u.a[w>>1].end&1);
    //         fprintf(stderr, "v_uId>>1: %u, v_uId&1: %u\n", v_uId>>1, v_uId&1);
    //         fprintf(stderr, "w_uId>>1: %u, w_uId&1: %u\n", w_uId>>1, w_uId&1);
    //         v = new_rtg_edges.a[k].ul>>32;
    //         w = new_rtg_edges.a[k].v;
    //         get_R_to_U(ruIndex, v>>1, &uId, &is_Unitig);
    //         fprintf(stderr, "v(read)>>1: %u, v(read)&1: %u, uId: %u, is_Unitig: %u\n", 
    //                                     v>>1, v&1, uId, is_Unitig);
    //         get_R_to_U(ruIndex, w>>1, &uId, &is_Unitig);
    //         fprintf(stderr, "w(read)>>1: %u, w(read)&1: %u, uId: %u, is_Unitig: %u\n", 
    //                                     w>>1, w&1, uId, is_Unitig);

            
    //     }
    // }
    /*********************************for debug**************************************/

    for (k = 0; k < new_utg_edges.n; k++)
    {
        v = new_utg_edges.a[k].ul>>32;
        w = new_utg_edges.a[k].v;

        if(nsg->seq[v>>1].del) continue;
        if(nsg->seq[w>>1].del) continue;
        ///if this edge is at a bubble
        if(nsg->seq_vis[v]!=0 && nsg->seq_vis[w^1]!=0)
        {
            ///fprintf(stderr, "v>>1: %u, w>>1: %u\n", v>>1, w>>1);
            continue;
        } 

        asg_arc_del(nsg, v, w, 1);
        asg_arc_del(nsg, w^1, v^1, 1);
        new_rtg_edges.a[k].del = 1;
    }
    asg_cleanup(nsg);
    free(nsg->seq_vis); nsg->seq_vis = NULL;

    /*********************************for debug**************************************/
    asg_arc_t* p = NULL;
    uint32_t e_occ = 0;
    for (k = 0; k < new_rtg_edges.n; k++)
    {
        if(new_rtg_edges.a[k].del) continue;
        p = asg_arc_pushp(r_g);
        *p = new_rtg_edges.a[k];
        e_occ++;
        if(keep_edges) kv_push(asg_arc_t, keep_edges->a, new_rtg_edges.a[k]);
    }

    if(e_occ != 0)
    {
        free(r_g->idx);
        r_g->idx = 0;
        r_g->is_srt = 0;
        asg_cleanup(r_g);
    }
    /*********************************for debug**************************************/

    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }

    if(i_ug == NULL) ma_ug_destroy(ug);
    kv_destroy(new_utg_edges);
    kv_destroy(new_rtg_edges);
    kv_destroy(rbub_edges);
}




///find contained read with longest overlap
ma_hit_t* get_best_no_coverage_read(ma_ug_t *ug, asg_t *r_g, ma_hit_t_alloc* reverse_sources, 
ma_sub_t *coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, uint32_t query,
uint32_t init_contain_uId)
{
    uint32_t qn = query>>1;
    int32_t r;
    asg_arc_t t;
    ma_hit_t *h = NULL, *return_h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t_alloc* x = &(reverse_sources[qn]);
    uint32_t i, is_Unitig, contain_uId;
    uint32_t maxOlen = 0;
    asg_t* nsg = ug->g;


    //scan all edges of qn
    for (i = 0; i < x->length; i++)
    {
        h = &(x->buffer[i]);
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);

        ///all of them cannot be removed
        if(sq->del || st->del || h->del) continue;
        if(r_g->seq[Get_qn(*h)].del || r_g->seq[Get_tn(*h)].del) continue;

        get_R_to_U(ruIndex, Get_tn(*h), &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) continue;
        if(nsg->seq[contain_uId].del) continue;
        ///contain_uId must be a unitig

        if(init_contain_uId != contain_uId) continue;

        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t);
        
        ///if it is a contained read, skip
        if(r < 0) continue;

        if((t.ul>>32) != query) continue;

        if(t.ol > maxOlen)
        {
            maxOlen = t.ol;
            return_h = h;
        }
    }

    return return_h;
}





int get_no_coverage_reads_chain_simple(ma_hit_t *h, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, R_to_U* ruIndex, ma_ug_t *ug, 
asg_t *r_g, int max_hang, int min_ovlp, uint32_t endRid, uint32_t uId, 
kvec_t_u32_warp* chain_buffer, kvec_asg_arc_t_warp* chain_edges, uint32_t* return_ava_cur, 
uint32_t* return_ava_ol, uint32_t* return_chainLen, uint32_t thresLen)
{
    (*return_chainLen) = (*return_ava_cur) = (*return_ava_ol) = (uint32_t)-1;
    uint32_t init_contain_uId = (uint32_t)-1, contain_uId, is_Unitig;
    uint32_t chainLen = 0, ava_cur, test_oLen;
    int ql, tl;
    int32_t r;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    asg_arc_t t;
    asg_t* nsg = ug->g;
    chain_buffer->a.n = 0;
    kv_push(uint32_t, chain_buffer->a, uId);
    if(chain_edges) chain_edges->a.n = 0;

    if(!h) return 0;
    sq = &(coverage_cut[Get_qn(*h)]);
    st = &(coverage_cut[Get_tn(*h)]);
    ///sq and st cannot be deleted
    ///for h, deleted or not does not matter
    if(sq->del || st->del) return 0;
    if(r_g->seq[Get_qn(*h)].del || r_g->seq[Get_tn(*h)].del) return 0;
    get_R_to_U(ruIndex, Get_tn(*h), &contain_uId, &is_Unitig);
    if(contain_uId == (uint32_t)-1 || is_Unitig != 1) return 0;
    if(nsg->seq[contain_uId].del) return 0;
    ql = sq->e - sq->s; tl = st->e - st->s;
    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
    if(r < 0) return 0;
    if((t.ul>>32) != endRid) return 0;
    if(chain_edges) kv_push(asg_arc_t, chain_edges->a, t);
    chainLen++;
    kv_push(uint32_t, chain_buffer->a, contain_uId);
    ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, reverse_sources, coverage_cut, 
            ruIndex, max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);
    //means find an aim
    if(ava_cur > 0)
    {
        (*return_ava_cur) = ava_cur;
        (*return_ava_ol) = test_oLen;
        (*return_chainLen) = chainLen;
        h = NULL;
        return 1;
    }
    else
    {
        return 0;
    }
    

    ///continue
    ///need to update h, endRid, chainLen, chain_buffer and chain_edges
    while (h)
    {
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);

        ///sq, st, and h cannot be deleted
        if(sq->del || st->del || h->del) break;
        if(r_g->seq[Get_qn(*h)].del || r_g->seq[Get_tn(*h)].del) break;

        get_R_to_U(ruIndex, Get_tn(*h), &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) break;
        if(nsg->seq[contain_uId].del) break;
        ///contain_uId must be a unitig

        if(init_contain_uId != (uint32_t)-1)
        {
            init_contain_uId = contain_uId;
        }
        else
        {
            if(init_contain_uId != contain_uId) break;
        }
         

        ql = sq->e - sq->s; tl = st->e - st->s;
        r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

        ///if st is contained in sq, or vice verse, skip
        if(r < 0) break;

        ///if sq and v are not in the same direction, skip
        ///endRid is (t.ul>>32), and t.v is a contained read
        if((t.ul>>32) != endRid) break;

        if(chain_edges) kv_push(asg_arc_t, chain_edges->a, t);
        chainLen++;
        kv_push(uint32_t, chain_buffer->a, contain_uId);
        ///endRid is (t.ul>>32), and t.v is a contained read
        ///find edges from t.v to existing unitigs
        ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, reverse_sources, coverage_cut, 
            ruIndex, max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);
        //means find an aim
        if(ava_cur > 0)
        {
            /**
            //do nothing if the chainLen == 1
            if(chainLen > 1)
            {
                asg_arc_t t_max;
                ma_hit_t_alloc* x = &(sources[(endRid>>1)]);
                uint32_t ava_ol_max = 0, ava_max = 0, k;
                for (k = 0; k < x->length; k++)
                {
                    h = &(x->buffer[k]);
                    sq = &(coverage_cut[Get_qn(*h)]);
                    st = &(coverage_cut[Get_tn(*h)]);
                    trio_flag = R_INF.trio_flag[Get_qn(*h)];

                    ///don't want to edges between different haps
                    non_trio_flag = (uint32_t)-1;
                    if(trio_flag == FATHER) non_trio_flag = MOTHER;
                    if(trio_flag == MOTHER) non_trio_flag = FATHER;
                    if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;
                    
                    ///just need deleted edges
                    ///sq might be deleted or not
                    if(!st->del) continue;
                    if(!h->del) continue;

                    ///tn must be contained in another existing read
                    get_R_to_U(ruIndex, Get_tn(*h), &contain_rId, &is_Unitig);
                    if(contain_rId == (uint32_t)-1 || is_Unitig == 1) continue;
                    if(r_g->seq[contain_rId].del) continue;

                    get_R_to_U(ruIndex, contain_rId, &contain_uId, &is_Unitig);
                    if(contain_uId == (uint32_t)-1 || is_Unitig != 1) continue;
                    if(nsg->seq[contain_uId].del) continue;
                    ///contain_uId must be a unitig

                    ql = sq->e - sq->s; tl = st->e - st->s;
                    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

                    ///if st is contained in sq, or vice verse, skip
                    if(r < 0) continue;

                    ///if sq and v are not in the same direction, skip
                    if((t.ul>>32) != endRid) continue;

                    ///pop the last contain_uId
                    chain_buffer->a.n--;
                    kv_push(uint32_t, chain_buffer->a, contain_uId);

                    ///note: t.v is a contained read
                    ///we need to find existing reads linked with w
                    ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);

                    if((ava_cur > ava_max) || (ava_cur == ava_max && test_oLen > ava_ol_max))
                    {
                        ava_max = ava_cur;
                        ava_ol_max = test_oLen;
                        t_max = t;
                    }
                }

                if(ava_max > 0)
                {
                    ava_cur = ava_max;
                    test_oLen = ava_ol_max;
                    if(chain_edges)
                    {
                        //pop the last edge
                        chain_edges->a.n--;
                        kv_push(asg_arc_t, chain_edges->a, t_max);
                    } 
                }
                else
                {
                    break;
                }
            }
            **/

            (*return_ava_cur) = ava_cur;
            (*return_ava_ol) = test_oLen;
            (*return_chainLen) = chainLen;
            h = NULL;
            return 1;
        }

        if(chainLen >= thresLen) break;

        ///endRid is (t.ul>>32), and t.v is a contained read
        ///haven't found a existing unitig from t.v
        ///check if t.v can link to a new contained read
        endRid = t.v; 
        h = get_best_no_coverage_read(ug, r_g, reverse_sources, coverage_cut, ruIndex, 
            max_hang, min_ovlp, endRid, init_contain_uId);
    }

    return 0;    
}


int get_no_coverage_reads_chain(ma_hit_t *h, ma_hit_t_alloc* sources, 
ma_hit_t_alloc* reverse_sources, ma_sub_t *coverage_cut, R_to_U* ruIndex, ma_ug_t *ug, 
asg_t *r_g, int max_hang, int min_ovlp, uint32_t endRid, uint32_t uId, 
kvec_t_u32_warp* chain_buffer, kvec_asg_arc_t_warp* chain_edges, uint32_t* return_ava_cur, 
uint32_t* return_ava_ol, uint32_t* return_chainLen, uint32_t thresLen)
{
    (*return_chainLen) = (*return_ava_cur) = (*return_ava_ol) = (uint32_t)-1;
    uint32_t init_contain_uId = (uint32_t)-1, contain_uId, is_Unitig;
    uint32_t chainLen = 0, ava_cur, test_oLen;
    int ql, tl;
    int32_t r;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    asg_arc_t t;
    asg_t* nsg = ug->g;
    chain_buffer->a.n = 0;
    kv_push(uint32_t, chain_buffer->a, uId);
    if(chain_edges) chain_edges->a.n = 0;

    if(!h) return 0;

    while(h)
    {
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);
        ///sq and st cannot be deleted
        ///for h, deleted or not does not matter
        if(sq->del || st->del) return 0;
        if(r_g->seq[Get_qn(*h)].del || r_g->seq[Get_tn(*h)].del) return 0;
        get_R_to_U(ruIndex, Get_tn(*h), &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) return 0;
        if(nsg->seq[contain_uId].del) return 0;
        ql = sq->e - sq->s; tl = st->e - st->s;
        r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);
        if(r < 0) return 0;
        if((t.ul>>32) != endRid) return 0;
        if(chain_edges) kv_push(asg_arc_t, chain_edges->a, t);
        chainLen++;
        kv_push(uint32_t, chain_buffer->a, contain_uId);
        ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, reverse_sources, coverage_cut, 
                ruIndex, max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);
        //means find an aim
        if(ava_cur > 0)
        {
            (*return_ava_cur) = ava_cur;
            (*return_ava_ol) = test_oLen;
            (*return_chainLen) = chainLen;
            h = NULL;
            return 1;
        }
        

        if(chainLen >= thresLen) break;

        ///endRid is (t.ul>>32), and t.v is a contained read
        ///haven't found a existing unitig from t.v
        ///check if t.v can link to a new contained read
        endRid = t.v; 
        h = get_best_no_coverage_read(ug, r_g, sources, coverage_cut, ruIndex, 
            max_hang, min_ovlp, endRid, init_contain_uId);
    }
    
    

    ///continue
    ///need to update h, endRid, chainLen, chain_buffer and chain_edges
    while (h)
    {
        sq = &(coverage_cut[Get_qn(*h)]);
        st = &(coverage_cut[Get_tn(*h)]);

        ///sq, st, and h cannot be deleted
        if(sq->del || st->del || h->del) break;
        if(r_g->seq[Get_qn(*h)].del || r_g->seq[Get_tn(*h)].del) break;

        get_R_to_U(ruIndex, Get_tn(*h), &contain_uId, &is_Unitig);
        if(contain_uId == (uint32_t)-1 || is_Unitig != 1) break;
        if(nsg->seq[contain_uId].del) break;
        ///contain_uId must be a unitig

        if(init_contain_uId != (uint32_t)-1)
        {
            init_contain_uId = contain_uId;
        }
        else
        {
            if(init_contain_uId != contain_uId) break;
        }
         

        ql = sq->e - sq->s; tl = st->e - st->s;
        r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

        ///if st is contained in sq, or vice verse, skip
        if(r < 0) break;

        ///if sq and v are not in the same direction, skip
        ///endRid is (t.ul>>32), and t.v is a contained read
        if((t.ul>>32) != endRid) break;

        if(chain_edges) kv_push(asg_arc_t, chain_edges->a, t);
        chainLen++;
        kv_push(uint32_t, chain_buffer->a, contain_uId);
        ///endRid is (t.ul>>32), and t.v is a contained read
        ///find edges from t.v to existing unitigs
        ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, reverse_sources, coverage_cut, 
            ruIndex, max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);
        //means find an aim
        if(ava_cur > 0)
        {
            /**
            //do nothing if the chainLen == 1
            if(chainLen > 1)
            {
                asg_arc_t t_max;
                ma_hit_t_alloc* x = &(sources[(endRid>>1)]);
                uint32_t ava_ol_max = 0, ava_max = 0, k;
                for (k = 0; k < x->length; k++)
                {
                    h = &(x->buffer[k]);
                    sq = &(coverage_cut[Get_qn(*h)]);
                    st = &(coverage_cut[Get_tn(*h)]);
                    trio_flag = R_INF.trio_flag[Get_qn(*h)];

                    ///don't want to edges between different haps
                    non_trio_flag = (uint32_t)-1;
                    if(trio_flag == FATHER) non_trio_flag = MOTHER;
                    if(trio_flag == MOTHER) non_trio_flag = FATHER;
                    if(R_INF.trio_flag[Get_tn(*h)] == non_trio_flag) continue;
                    
                    ///just need deleted edges
                    ///sq might be deleted or not
                    if(!st->del) continue;
                    if(!h->del) continue;

                    ///tn must be contained in another existing read
                    get_R_to_U(ruIndex, Get_tn(*h), &contain_rId, &is_Unitig);
                    if(contain_rId == (uint32_t)-1 || is_Unitig == 1) continue;
                    if(r_g->seq[contain_rId].del) continue;

                    get_R_to_U(ruIndex, contain_rId, &contain_uId, &is_Unitig);
                    if(contain_uId == (uint32_t)-1 || is_Unitig != 1) continue;
                    if(nsg->seq[contain_uId].del) continue;
                    ///contain_uId must be a unitig

                    ql = sq->e - sq->s; tl = st->e - st->s;
                    r = ma_hit2arc(h, ql, tl, max_hang, asm_opt.max_hang_rate, min_ovlp, &t);

                    ///if st is contained in sq, or vice verse, skip
                    if(r < 0) continue;

                    ///if sq and v are not in the same direction, skip
                    if((t.ul>>32) != endRid) continue;

                    ///pop the last contain_uId
                    chain_buffer->a.n--;
                    kv_push(uint32_t, chain_buffer->a, contain_uId);

                    ///note: t.v is a contained read
                    ///we need to find existing reads linked with w
                    ava_cur = get_num_edges2existing_nodes_advance(ug, r_g, sources, coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t.v, &test_oLen, chain_buffer->a.a, chain_buffer->a.n, 1);

                    if((ava_cur > ava_max) || (ava_cur == ava_max && test_oLen > ava_ol_max))
                    {
                        ava_max = ava_cur;
                        ava_ol_max = test_oLen;
                        t_max = t;
                    }
                }

                if(ava_max > 0)
                {
                    ava_cur = ava_max;
                    test_oLen = ava_ol_max;
                    if(chain_edges)
                    {
                        //pop the last edge
                        chain_edges->a.n--;
                        kv_push(asg_arc_t, chain_edges->a, t_max);
                    } 
                }
                else
                {
                    break;
                }
            }
            **/

            (*return_ava_cur) = ava_cur;
            (*return_ava_ol) = test_oLen;
            (*return_chainLen) = chainLen;
            h = NULL;
            return 1;
        }

        if(chainLen >= thresLen) break;

        ///endRid is (t.ul>>32), and t.v is a contained read
        ///haven't found a existing unitig from t.v
        ///check if t.v can link to a new contained read
        endRid = t.v; 
        h = get_best_no_coverage_read(ug, r_g, reverse_sources, coverage_cut, ruIndex, 
            max_hang, min_ovlp, endRid, init_contain_uId);
    }

    return 0;    
}


uint32_t check_if_no_coverage(asg_t *r_g, asg_arc_t* list, uint32_t list_n)
{
    if(list_n < 2) return 0;
    uint32_t i, v, w, u_n = list_n - 1, nv, l, k;
    long long totalLen = 0;
    asg_arc_t *av = NULL;
    for (i = 0; i < u_n - 1; i++)
    {
        v = list[i].v;
        w = list[i+1].v; 
        av = asg_arc_a(r_g, v);
        nv = asg_arc_n(r_g, v);
        l = 0;
        
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == w) 
            {
                l = asg_arc_len(av[k]);
                break;
            }
        }

        if(k == nv) fprintf(stderr, "####ERROR\n");
        totalLen += l;
    }

    if(i < u_n)
    {
        v = list[i].v;
        l = r_g->seq[v>>1].len;
        totalLen += l;
    }

    if(totalLen > list[0].ol + list[list_n - 1].ol) return 1;

    return 0;
}

uint32_t create_fake_read(asg_t *r_g, asg_arc_t* list, uint32_t list_n)
{
    if(list_n < 2) return 0;
    ma_utg_t* u = NULL;
    uint32_t v, w, l, i, k, nv, totalLen = 0;
    asg_arc_t *av = NULL;
    u = asg_F_seq_set(r_g, r_g->n_seq);
    if(u == NULL) return 0;
    u->n = u->m = list_n - 1;
    if(u->a) free(u->a);
    u->a = (uint64_t*)malloc(sizeof(uint64_t)*u->n);

    /*****************for debug**********************/
    // fprintf(stderr, "list_n: %u\n", list_n);
    // for (i = 0; i < list_n; i++)
    // {
    //     fprintf(stderr, "(%u) v: %u, dir: %u, w: %u, w&1: %u, len: %u, ol: %u\n", 
    //             i, (uint32_t)(list[i].ul>>33), (uint32_t)((list[i].ul>>32)&1), 
    //             list[i].v>>1, list[i].v&1, (uint32_t)list[i].ul, list[i].ol);
    // }
    /*****************for debug**********************/


    for (i = 0; i < u->n; i++)
    {
        u->a[i] = list[i].v;
        u->a[i] = u->a[i] << 32;
    }

    /*****************for debug**********************/
    // fprintf(stderr, "u->n: %u\n", u->n);
    // for (i = 0; i < u->n; i++)
    // {
    //     fprintf(stderr, "(%u) v: %u, dir: %u\n", 
    //                 i, (uint32_t)(u->a[i]>>33), (uint32_t)((u->a[i]>>32)&1));
    // }
    /*****************for debug**********************/
    
    
    for (i = 0; i < u->n - 1; i++)
    {
        v = (uint64_t)(u->a[i])>>32;
        w = (uint64_t)(u->a[i + 1])>>32; 
        av = asg_arc_a(r_g, v);
        nv = asg_arc_n(r_g, v);
        l = 0;
        
        for (k = 0; k < nv; k++)
        {
            if(av[k].del) continue;
            if(av[k].v == w) 
            {
                l = asg_arc_len(av[k]);
                break;
            }
        }

        if(k == nv) fprintf(stderr, "####ERROR\n");
        u->a[i] = v; u->a[i] = u->a[i]<<32; u->a[i] = u->a[i] | (uint64_t)(l);
        totalLen += l;
    }

    if(i < u->n)
    {
        v = (uint64_t)(u->a[i])>>32;
        l = r_g->seq[v>>1].len;
        u->a[i] = v;
        u->a[i] = u->a[i]<<32;
        u->a[i] = u->a[i] | (uint64_t)(l);
        totalLen += l;
    }

    /*****************for debug**********************/
    // fprintf(stderr, "u->n: %u\n", u->n);
    // for (i = 0; i < u->n; i++)
    // {
    //     fprintf(stderr, "(%u) v: %u, dir: %u, len: %u\n", 
    //                 i, (uint32_t)(u->a[i]>>33), (uint32_t)((u->a[i]>>32)&1), (uint32_t)(u->a[i]));
    // }
    /*****************for debug**********************/

    

    u->len = totalLen;
    u->circ = 0;
    u->start = list[0].ol;
    u->end = u->len - list[list_n - 1].ol;
    totalLen = u->len - list[0].ol - list[list_n - 1].ol;
    ///u->len = totalLen;

    /*****************for debug**********************/
    // fprintf(stderr, "u->len: %u, u->start: %u, u->end: %u\n", u->len, u->start, u->end);
    // fprintf(stderr, "totalLen: %u, list[0].ol: %u, list[list_n - 1].ol: %u\n", 
    // totalLen, list[0].ol, list[list_n - 1].ol);
    /*****************for debug**********************/

    /*****************for debug**********************/
    /**
    (*coverage_cut) = (ma_sub_t*)realloc((*coverage_cut), (r_g->n_seq+1)*sizeof(ma_sub_t));
    (*coverage_cut)[r_g->n_seq].del = 0;
    (*coverage_cut)[r_g->n_seq].c = PRIMARY_LABLE;
    (*coverage_cut)[r_g->n_seq].s = 0;
    (*coverage_cut)[r_g->n_seq].e = totalLen;
    **/
    asg_seq_set(r_g, r_g->n_seq, totalLen, 0);
    /*****************for debug**********************/
    return 1;     
}


inline void generate_edge(asg_t *r_g, uint32_t v, uint32_t w, uint32_t ol, uint32_t del,
uint32_t strong, uint32_t el, uint32_t no_l_indel, asg_arc_t* t)
{
    uint64_t l = r_g->seq[v>>1].len - ol;
    t->ul = v; t->ul = t->ul << 32; t->ul = t->ul|l;
    t->v = w;
    t->ol = ol;
    t->del = del;
    t->strong = strong;
    t->el = el;
    t->no_l_indel = no_l_indel;
}
///chainLenThres is used to avoid circle
void rescue_no_coverage_aggressive(asg_t *r_g, ma_hit_t_alloc* sources_count, 
ma_hit_t_alloc* reverse_source, ma_sub_t **coverage_cut, R_to_U* ruIndex, int max_hang, int min_ovlp, 
long long bubble_dist, uint32_t chainLenThres)
{
    uint32_t n_vtx, v, k, kv, is_Unitig, uId, rId, endRid, w;
    uint32_t ava_max, ava_ol_max, ava_min_chain, ava_cur, ava_chainLen, test_oLen, is_update;
    uint64_t interval_beg, intervalLen;
    asg_t* nsg = NULL;
    ma_utg_t* nsu = NULL;
    asg_arc_t t, t_max, r_edge;
    t_max.v = t_max.ul = t.v = t.ul = (uint32_t)-1;
    ma_hit_t *h = NULL, *h_max = NULL;
    ma_hit_t_alloc* x = NULL;
    ma_ug_t* ug = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    int32_t r;
    
    kvec_t(asg_arc_t) new_utg_edges;
    kv_init(new_utg_edges);

    kvec_t(asg_arc_t) new_rtg_edges;
    kv_init(new_rtg_edges);

    kvec_t(asg_arc_t) hap_edges;
    kv_init(hap_edges);

    kvec_t(asg_arc_t) rbub_edges;
    kv_init(rbub_edges);

    kvec_t_u64_warp u_vecs;
    kv_init(u_vecs.a);

    kvec_t_u64_warp intervals;
    kv_init(intervals.a);

    kvec_t_u32_warp chain_buffer;
    kv_init(chain_buffer.a);

    kvec_asg_arc_t_warp chain_edges;
    kv_init(chain_edges.a);


    ug = ma_ug_gen(r_g);
    nsg = ug->g;
    
    for (v = 0; v < nsg->n_seq; v++)
    {
        uId = v;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsg->seq[v].del) continue;
        nsg->seq[v].c = PRIMARY_LABLE;
        for (k = 0; k < nsu->n; k++)
        {
            rId = nsu->a[k]>>33;
            set_R_to_U(ruIndex, rId, uId, 1);
        }
    }




    n_vtx = nsg->n_seq * 2;
    for (v = 0; v < n_vtx; v++)
    {
        uId = v>>1;
        if(nsg->seq[uId].del) continue;
        nsu = &(ug->u.a[uId]);
        if(nsu->m == 0) continue;
        if(nsu->circ || nsu->start == UINT32_MAX || nsu->end == UINT32_MAX) continue;
        if(get_real_length(nsg, v, NULL) != 0) continue;
        ///we probably don't need this line
        ///if(get_real_length(nsg, v^1, NULL) == 0) continue;
        if(v&1)
        {
            endRid = nsu->start^1;
        }
        else
        {
            endRid = nsu->end^1;
        }


        x = &(sources_count[(endRid>>1)]);
        for (k = 0; k < x->length; k++)
        {
            ///h is the edge of endRid
            h = &(x->buffer[k]);
            sq = &((*coverage_cut)[Get_qn(*h)]);
            st = &((*coverage_cut)[Get_tn(*h)]);
            ///sq, st, and h cannot be deleted
            if(sq->del || st->del || h->del) continue;
            r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                            asm_opt.max_hang_rate, min_ovlp, &t);
            if(r < 0) continue;
            if((t.ul>>32) != endRid) continue;
            break;
        }
        ///means there is >0 edges
        if(k != x->length) continue;

        ///x is the end read of a tip
        ///find all overlap of x
        x = &(reverse_source[(endRid>>1)]);
        ava_ol_max = ava_max = 0; ava_min_chain = (uint32_t)-1;
        h_max = NULL;
        for (k = 0; k < x->length; k++)
        {
            h = &(x->buffer[k]);
            ///means we found a contained read
            if(get_no_coverage_reads_chain_simple(h, sources_count, reverse_source, *coverage_cut, 
            ruIndex, ug, r_g, max_hang, min_ovlp, endRid, uId, &chain_buffer, NULL, &ava_cur, 
            &test_oLen, &ava_chainLen, chainLenThres))
            {
                
                is_update = 0;
                if(ava_chainLen < ava_min_chain)
                {
                    is_update = 1;
                }
                else if(ava_chainLen == ava_min_chain)
                {
                    if(ava_cur > ava_max)
                    {
                        is_update = 1;
                    }
                    else if(ava_cur == ava_max && test_oLen > ava_ol_max)
                    {
                        is_update = 1;
                    }
                }
                if(is_update)
                {
                    ava_min_chain = ava_chainLen;
                    ava_max = ava_cur;
                    ava_ol_max = test_oLen;
                    h_max = h;
                }
            }
        }





        if(ava_max > 0)
        {

            //get all edges between contained reads
            get_no_coverage_reads_chain_simple(h_max, sources_count, reverse_source, 
            *coverage_cut, ruIndex, ug, r_g, max_hang, min_ovlp, endRid, uId, &chain_buffer, 
            &chain_edges, &ava_cur, &test_oLen, &ava_chainLen, chainLenThres);
            if(chain_edges.a.n < 1) continue;
            ///the last cantained read
            t_max = chain_edges.a.a[chain_edges.a.n-1];

            k = 0; rbub_edges.n = 0;
            ///edges from the last contained read to other unitigs
            while(get_edge2existing_node_advance(ug, r_g, reverse_source, *coverage_cut, ruIndex, 
                    max_hang, min_ovlp, t_max.v, chain_buffer.a.a, chain_buffer.a.n, &k, &r_edge, 1))
            {
               kv_push(asg_arc_t, rbub_edges, r_edge);
            }

            ///need to do transitive reduction
            ///note here is different to standard transitive reduction
            minor_transitive_reduction(ug, ruIndex, rbub_edges.a, rbub_edges.n);

            for (k = 0, kv = 0; k < rbub_edges.n; k++)
            {
                t = rbub_edges.a[k];
                if(t.del) continue;
                kv++;
            }

            if(kv != 1) continue;

            interval_beg = hap_edges.n;
            
            //just saves all nodes here
            ///note that only the first edge is generated from reverse_source, 
            ///others are generated from source
            for (k = 0; k < chain_edges.a.n; k++)
            {
                kv_push(asg_arc_t, hap_edges, chain_edges.a.a[k]);
            }

            for (k = 0; k < rbub_edges.n; k++)
            {
                t = rbub_edges.a[k];
                if(t.del) continue;
                kv_push(asg_arc_t, hap_edges, t);
            }


            if(check_if_no_coverage(r_g, hap_edges.a+interval_beg, hap_edges.n - interval_beg) == 0)
            {
                hap_edges.n = interval_beg;
                continue;
            }

            ///there is just one edge
            ///connect multiple edges is too difficult
            for (k = 0; k < rbub_edges.n; k++)
            {
                t = rbub_edges.a[k];
                if(t.del) continue;

                w = get_corresponding_uId(ug, ruIndex, t.v);
                if(w == ((uint32_t)(-1))) continue;
                kv_push(asg_arc_t, new_rtg_edges, t);

                get_edge_from_source(reverse_source, *coverage_cut, ruIndex, max_hang, min_ovlp, 
                (t.v^1), ((t.ul>>32)^1), &t);
                kv_push(asg_arc_t, new_rtg_edges, t);


                ///for unitig graph
                asg_append_edges_to_srt(nsg, v, ug->u.a[v>>1].len, w, 0, 0, 0, 0);
                asg_append_edges_to_srt(nsg, w^1, ug->u.a[w>>1].len, v^1, 0, 0, 0, 0);

                t.ul = v; t.ul = t.ul<<32; t.v = w;
                kv_push(asg_arc_t, new_utg_edges, t);

                t.ul = w^1; t.ul = t.ul<<32; t.v = v^1;
                kv_push(asg_arc_t, new_utg_edges, t);

            }

            kv_push(uint64_t, u_vecs.a, v);

            intervalLen = hap_edges.n - interval_beg;
            interval_beg = interval_beg << 32;
            interval_beg = interval_beg | intervalLen;
            kv_push(uint64_t, intervals.a, interval_beg);
        }
    }

    nsg->seq_vis = (uint8_t*)calloc(nsg->n_seq*2, sizeof(uint8_t));    
    lable_all_bubbles(nsg, bubble_dist);


    for (k = 0; k < new_utg_edges.n; k++)
    {
        v = new_utg_edges.a[k].ul>>32;
        w = new_utg_edges.a[k].v;

        if(nsg->seq[v>>1].del) continue;
        if(nsg->seq[w>>1].del) continue;
        ///if this edge is at a bubble
        if(nsg->seq_vis[v]!=0 && nsg->seq_vis[w^1]!=0)
        {
            continue;
        } 

        asg_arc_del(nsg, v, w, 1);
        asg_arc_del(nsg, w^1, v^1, 1);
        new_rtg_edges.a[k].del = 1;
    }
    free(nsg->seq_vis); nsg->seq_vis = NULL;


    new_utg_edges.n = 0;
    asg_arc_t* list = NULL;
    uint32_t nextNode, curNode, pre_num_nodes = r_g->n_seq;
    ///u_vecs saves the unitig Id
    for (k = 0; k < u_vecs.a.n; k++)
    {
        w = (uint32_t)u_vecs.a.a[k];
        ///do nothing
        if(get_real_length(r_g, w, NULL) == 0) continue;

        interval_beg = intervals.a.a[k]>>32;
        intervalLen = intervals.a.a[k] & ((uint64_t)0xffffffff);
        list = hap_edges.a + interval_beg;
        if(intervalLen < 2) continue;///fprintf(stderr, "No enough edges\n");
        fprintf(stderr, "\n(%u) w>>1: %u, w&1: %u\n", k, w>>1, w&1);
        ///continue;


        if(create_fake_read(r_g, list, intervalLen) == 0) continue;


        curNode = list[0].ul>>32;
        nextNode = (r_g->n_seq - 1)<<1;
        /*****************for debug**********************/
        generate_edge(r_g, curNode, nextNode, 0, 0, 0, 0, 0, &t);
        kv_push(asg_arc_t, new_utg_edges, t);

        generate_edge(r_g, nextNode^1, curNode^1, 0, 0, 0, 0, 0, &t);
        kv_push(asg_arc_t, new_utg_edges, t);
        /*****************for debug**********************/
        // fprintf(stderr, "curNode>>1: %u, curNode&1: %u, nextNode>>1: %u, nextNode&1: %u\n", 
        // curNode>>1, curNode&1, nextNode>>1, nextNode&1);


        curNode = (r_g->n_seq - 1)<<1;
        nextNode = list[intervalLen - 1].v;
        /*****************for debug**********************/
        generate_edge(r_g, curNode, nextNode, 0, 0, 0, 0, 0, &t);
        kv_push(asg_arc_t, new_utg_edges, t);

        generate_edge(r_g, nextNode^1, curNode^1, 0, 0, 0, 0, 0, &t);
        kv_push(asg_arc_t, new_utg_edges, t);
        /*****************for debug**********************/
        // fprintf(stderr, "curNode>>1: %u, curNode&1: %u, nextNode>>1: %u, nextNode&1: %u\n", 
        // curNode>>1, curNode&1, nextNode>>1, nextNode&1);
    }

    


    asg_arc_t* p = NULL;
    for (k = 0; k < new_utg_edges.n; k++)
    {
        p = asg_arc_pushp(r_g);
        *p = new_utg_edges.a[k];
    }

    
    free(r_g->idx);
    r_g->idx = 0;
    r_g->is_srt = 0;
    asg_cleanup(r_g);

    if(r_g->n_seq > pre_num_nodes)
    {
        (*coverage_cut) = (ma_sub_t*)realloc((*coverage_cut), r_g->n_seq*sizeof(ma_sub_t));
        R_INF.trio_flag = (uint8_t*)realloc(R_INF.trio_flag, r_g->n_seq*sizeof(uint8_t));
        ruIndex->index = (uint32_t*)realloc(ruIndex->index, r_g->n_seq*sizeof(uint32_t));
        ruIndex->len = r_g->n_seq;
        for (k = pre_num_nodes; k < r_g->n_seq; k++)
        {
            ruIndex->index[k] = (uint32_t)-1;
            R_INF.trio_flag[k] = AMBIGU;
            (*coverage_cut)[k].del = 0;
            (*coverage_cut)[k].c = PRIMARY_LABLE;
            (*coverage_cut)[k].s = 0;
            (*coverage_cut)[k].e = r_g->seq[k].len;
        }
    }
    

    for (v = 0; v < ruIndex->len; v++)
    {
        get_R_to_U(ruIndex, v, &uId, &is_Unitig);
        if(is_Unitig == 1) ruIndex->index[v] = (uint32_t)-1;
    }


    ma_ug_destroy(ug);
    kv_destroy(new_utg_edges);
    kv_destroy(hap_edges);
    kv_destroy(new_rtg_edges);
    kv_destroy(rbub_edges);
    kv_destroy(u_vecs.a);
    kv_destroy(chain_buffer.a);
    kv_destroy(chain_edges.a);
    kv_destroy(intervals.a);
}


void fix_binned_reads(ma_hit_t_alloc* paf, uint64_t n_read, ma_sub_t* coverage_cut)
{
    double startTime = Get_T();
    uint64_t i, binned_flag, binned_reads = 0, binned_error_reads = 0, reduce = 1;
    ma_sub_t max_left, max_right;
    while (reduce != 0)
    {    
        reduce = 0;
        binned_reads = 0;
        for (i = 0; i < n_read; ++i) 
        {
            if(coverage_cut[i].del) continue;
            binned_flag = R_INF.trio_flag[i];
            if(binned_flag == AMBIGU) continue;
            binned_reads++;
            max_left.s = max_right.s = Get_READ_LENGTH(R_INF,i);
            max_left.e = max_right.e = 0;
            collect_sides_trio(&(paf[i]), Get_READ_LENGTH(R_INF,i), &max_left, &max_right, binned_flag);
            collect_contain_trio(&(paf[i]), NULL, Get_READ_LENGTH(R_INF,i), &max_left, &max_right, 0.1,
            binned_flag);
            if(max_left.e > max_right.s)
            {
                R_INF.trio_flag[i] = AMBIGU;
                binned_error_reads++;
                reduce++;
                binned_reads--;
            } 
            /**
            if(max_left.e > max_right.s)
            {
                fprintf(stderr, "\ni: %lu, reference: %.*s, len: %lu, %c\n", 
                    i, (int)Get_NAME_LENGTH(R_INF, i), Get_NAME(R_INF, i), 
                    Get_READ_LENGTH(R_INF, i), "apmaaa"[binned_flag]);

                for (uint64_t j = 0; j < paf[i].length; j++)
                {
                    if(R_INF.trio_flag[Get_tn(paf[i].buffer[j])] == binned_flag)
                    {
                        fprintf(stderr, "###Compatible, ");
                    }
                    else
                    {
                        fprintf(stderr, "***Conflict, ");
                    }

                    fprintf(stderr, "%c, qs: %u, qe: %u, ts: %u, te: %u, len: %lu, %.*s\n", 
                    "apmaaa"[R_INF.trio_flag[Get_tn(paf[i].buffer[j])]],
                    Get_qs(paf[i].buffer[j]),
                    Get_qe(paf[i].buffer[j]),
                    Get_ts(paf[i].buffer[j]),
                    Get_te(paf[i].buffer[j]),
                    Get_READ_LENGTH(R_INF, Get_tn(paf[i].buffer[j])),
                    (int)Get_NAME_LENGTH(R_INF, Get_tn(paf[i].buffer[j])), 
                    Get_NAME(R_INF, Get_tn(paf[i].buffer[j])));
                }
            }
            **/
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    fprintf(stderr, "n_read: %lu, binned_reads: %lu, binned_error_reads: %lu\n", 
    (unsigned long)n_read, (unsigned long)binned_reads, (unsigned long)binned_error_reads);
}


void print_binned_reads(ma_hit_t_alloc* paf, uint64_t n_read, ma_sub_t* coverage_cut)
{
    uint64_t i, j, binned_flag;

    for (i = 0; i < n_read; ++i)
    {
        if(coverage_cut!=NULL && coverage_cut[i].del) continue;
        binned_flag = R_INF.trio_flag[i];
        fprintf(stdout, "\n***i: %u, binned_flag: %u\n", (uint32_t)i, (uint32_t)binned_flag);
        if(coverage_cut!=NULL)
        {
            fprintf(stdout, "coverage_cut[i].c: %u, coverage_cut[i].s: %u, coverage_cut[i].e: %u, coverage_cut[i].e: %u\n",
            coverage_cut[i].c, coverage_cut[i].s, coverage_cut[i].e, coverage_cut[i].del);
        }
        for (j = 0; j < paf[i].length; j++)
        {
            if(paf[i].buffer[j].del) continue;
            ///fprintf(stdout, "j: %u\n", (uint32_t)j);
            fprintf(stdout, "qn: %u, qs: %u, qe: %u, tn: %u, ts: %u, te: %u\n", 
            (uint32_t)Get_qn(paf[i].buffer[j]), (uint32_t)Get_qs(paf[i].buffer[j]), (uint32_t)Get_qe(paf[i].buffer[j]),
            (uint32_t)Get_tn(paf[i].buffer[j]), (uint32_t)Get_ts(paf[i].buffer[j]), (uint32_t)Get_te(paf[i].buffer[j]));
        }
    }
}



void debug_ma_hit_t(ma_hit_t_alloc* sources, ma_sub_t *coverage_cut, long long num_sources,
int max_hang, int min_ovlp)
{
    double startTime = Get_T();

    long long i, j, index;
    uint32_t qn, tn;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    ma_hit_t *h = NULL;
    int32_t r_0, r_1;
    asg_arc_t t_0, t_1;
    fprintf(stderr, "start!\n");

    for (i = 0; i < num_sources; i++)
    {

        for (j = 0; j < sources[i].length; j++)
        {
            qn = Get_qn(sources[i].buffer[j]);
            tn = Get_tn(sources[i].buffer[j]);


            ///if(sources[i].buffer[j].del) continue;

            index = get_specific_overlap(&(sources[tn]), tn, qn);
            if(index == -1)
            {
                fprintf(stderr, "[E::%s]ERROR 0: qn: %u, tn: %u\n", __func__, qn, tn );
                continue;
            } 



            if(sources[i].buffer[j].del != sources[tn].buffer[index].del)
            {
                fprintf(stderr, "[E::%s]ERROR 2: qn: %u, tn: %u\n", __func__, qn, tn);
            }

            h = &(sources[i].buffer[j]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            r_0 = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_0);


            h = &(sources[tn].buffer[index]);
            sq = &(coverage_cut[Get_qn(*h)]);
            st = &(coverage_cut[Get_tn(*h)]);
            r_1 = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, 
                                asm_opt.max_hang_rate, min_ovlp, &t_1);
            t_0.ul=t_0.v=t_1.ul=t_1.v = 0;
            if(r_0 < 0 && r_1 < 0) continue;
            if((t_0.ul>>32) != (t_1.v^1)) fprintf(stderr, "[E::%s]ERROR 3: qn: %u, tn: %u\n", __func__, qn, tn);
            if((t_1.ul>>32) != (t_0.v^1)) fprintf(stderr, "[E::%s]ERROR 4: qn: %u, tn: %u\n", __func__, qn, tn);
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] takes %0.2fs\n\n", __func__, Get_T()-startTime);
    }
}


#define mahit_tn_key(a) ((a).tn)
KRADIX_SORT_INIT(pafbuffer_sort, ma_hit_t, mahit_tn_key, member_size(ma_hit_t, tn))

typedef struct{
    long long n_read;
    ma_hit_t_alloc *x;

}pafbuffer_sort_tn_aux_t;
static void sort_paf_buffers_by_targetID_worker(void *data, long i_read, int tid){  // callback for kt_for
    pafbuffer_sort_tn_aux_t *aux = (pafbuffer_sort_tn_aux_t*)data;
    uint32_t buf_l = aux->x[i_read].length;

    if (buf_l<LINEAR_BF_BREAKEVEN_POINT) return;  // will use linear search
    if (aux->x[i_read].buffer_sorted_by_tn==(buf_l-1)) return;  // already sorted

    ma_hit_t *buf = aux->x[i_read].buffer;
    radix_sort_pafbuffer_sort(buf, buf+buf_l);
    aux->x[i_read].buffer_sorted_by_tn = buf_l-1;
    return;
}
void sort_paf_buffers_by_targetID(ma_hit_t_alloc *source, long long n_read){
    // FUNC
    //     Sort each x[qn].buffer by its x[qn].buffer[i].tn values.
    //     Used to enable hamt's binary searches.
    double startTime = Get_T();
    int n_cpu = asm_opt.thread_num;
    pafbuffer_sort_tn_aux_t aux;
    aux.x = source;
    aux.n_read = n_read;

    kt_for(n_cpu, sort_paf_buffers_by_targetID_worker, &aux, n_read);

    // fprintf(stderr, "[T::%s] used %.2f s\n", __func__, Get_T()-startTime);
}
void sort_one_paf_buffer_by_targetID(ma_hit_t_alloc *x){
    if (x->length<LINEAR_BF_BREAKEVEN_POINT) return;
    ma_hit_t *buf = x->buffer;
    radix_sort_pafbuffer_sort(x->buffer, x->buffer+x->length);
    x->buffer_sorted_by_tn = 1;
}
void reset_paf_buffers_sorting_marks_targetID(ma_hit_t_alloc* sources, long long n_read){
    for (long long i=0; i<n_read; i++){
        sources[i].buffer_sorted_by_tn = 0;
    }
}


#define mahit_qns_key(a) ((a).qns)
KRADIX_SORT_INIT(pafbuffer_sort_qns, ma_hit_t, mahit_qns_key, member_size(ma_hit_t, qns))
typedef struct{
    long long n_read;
    ma_hit_t_alloc *x;
}pafbuffer_sort_qns_aux_t;
static void sort_paf_buffers_by_qns_worker(void *data, long i_read, int tid){  // callback for kt_for
    pafbuffer_sort_qns_aux_t *aux = (pafbuffer_sort_qns_aux_t*)data;
    uint32_t buf_l = aux->x[i_read].length;
    ma_hit_t *buf = aux->x[i_read].buffer;
    radix_sort_pafbuffer_sort_qns(buf, buf+buf_l);
}
void sort_paf_buffers_by_qns(ma_hit_t_alloc *source, long long n_read){
    // FUNC
    //     Used to restore the ordering (key is qns) after binary search.
    double startTime = Get_T();
    int n_cpu = asm_opt.thread_num;
    pafbuffer_sort_qns_aux_t aux;
    aux.x = source;
    aux.n_read = n_read;
    kt_for(n_cpu, sort_paf_buffers_by_qns_worker, &aux, n_read);
}


void clean_graph(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, long long mini_overlap_length, 
long long max_hang_length, long long clean_round, long long gap_fuzz,
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, char* output_file_name, 
long long bubble_dist, int read_graph, R_to_U* ruIndex, asg_t **sg_ptr,   // the **sg_ptr is a NULL pointer
ma_sub_t **coverage_cut_ptr, int debug_g)
{
		ma_sub_t *coverage_cut = *coverage_cut_ptr;
	asg_t *sg = *sg_ptr;

    if(debug_g) goto debug_gfa;

    ///just for debug
    renew_graph_init(sources, reverse_sources, sg, coverage_cut, ruIndex, n_read);

    ///it's hard to say which function is better       
    ///normalize_ma_hit_t_single_side(sources, n_read);
    normalize_ma_hit_t_single_side_advance(sources, n_read);
    normalize_ma_hit_t_single_side_advance(reverse_sources, n_read);

    if (ha_opt_triobin(&asm_opt))
    {
        drop_edges_by_trio(sources, n_read);
    }
    else
    {
        memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads*sizeof(uint8_t));
    }
    
    ///print_binned_reads(sources, n_read, coverage_cut);

    clean_weak_ma_hit_t(sources, reverse_sources, n_read);
    ///ma_hit_sub is just use to init coverage_cut,
    ///it seems we do not need ma_hit_cut & ma_hit_flt
    ma_hit_sub(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    detect_chimeric_reads(sources, n_read, readLen, coverage_cut, asm_opt.max_ov_diff_final * 2.0);
    ma_hit_cut(sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ///print_binned_reads(sources, n_read, coverage_cut);
    ma_hit_flt(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);

    ///fix_binned_reads(sources, n_read, coverage_cut);
    ///just need to deal with trio here
    ma_hit_contained_advance(sources, n_read, coverage_cut, ruIndex, max_hang_length, mini_overlap_length);
    
    sg = ma_sg_gen(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);
    asg_arc_del_trans(sg, gap_fuzz);

    asm_opt.coverage = get_coverage(sources, coverage_cut, n_read);

    if(VERBOSE >= 1)
    {
        char* unlean_name = (char*)malloc(strlen(output_file_name)+25);
        sprintf(unlean_name, "%s.unclean", output_file_name);
        output_read_graph(sg, coverage_cut, unlean_name, n_read);
        free(unlean_name);
    }

    asg_cut_tip(sg, asm_opt.max_short_tip);
    ///debug_info_of_specfic_node("m64062_190803_042216/15205346/ccs", sg, "inner_1");
    ///drop_inexact_edegs_at_bubbles(sg, bubble_dist);
    
    if(clean_round > 0)
    {
        double cut_step;
        if(clean_round == 1)
        {
            cut_step = max_ovlp_drop_ratio;
        }
        else
        {
            cut_step = (max_ovlp_drop_ratio - min_ovlp_drop_ratio) / (clean_round - 1);
        }
        double drop_ratio = min_ovlp_drop_ratio;
        int i = 0;
        for (i = 0; i < clean_round; i++, drop_ratio += cut_step)
        {
            if(drop_ratio > max_ovlp_drop_ratio)
            {
                drop_ratio = max_ovlp_drop_ratio;
            }

            if(VERBOSE >= 1)
            {
                fprintf(stderr, "\n\n**********%d-th round drop: drop_ratio = %f**********\n", 
                i, drop_ratio);
            }

            ///just topological clean
            pre_clean(sources, coverage_cut, sg, bubble_dist);
            ///asg_arc_del_orthology(sg, reverse_sources, drop_ratio, asm_opt.max_short_tip);
            // asg_arc_del_orthology_multiple_way(sg, reverse_sources, drop_ratio, asm_opt.max_short_tip);
            // asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            //reomve edge between two chromesomes
            //this node must be a single read
            asg_arc_del_false_node(sg, sources, asm_opt.max_short_tip);
            asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/
            /****************************may have bugs********************************/

            ///asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_identify_simple_bubbles_multi(sg, 0);
            ///asg_arc_del_short_diploid_unclean_exact(sg, drop_ratio, sources);
            if (ha_opt_triobin(&asm_opt))
            {
                asg_arc_del_short_diploid_by_exact_trio(sg, asm_opt.max_short_tip, sources);
            }
            else
            {
                asg_arc_del_short_diploid_by_exact(sg, asm_opt.max_short_tip, sources);
            }
            asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            if (ha_opt_triobin(&asm_opt))
            {
                asg_arc_del_short_diploid_by_length_trio(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
                asm_opt.max_short_tip, 1, 1, 0, 0, ruIndex);
            }
            else
            {
                asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
                asm_opt.max_short_tip, 1, 1, 0, 0, ruIndex);
            }
            asg_cut_tip(sg, asm_opt.max_short_tip);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, 
            asm_opt.max_short_tip, ruIndex);     

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_complex_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, asm_opt.max_short_tip);

            asg_cut_tip(sg, asm_opt.max_short_tip);
        }
    }

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "\n\n**********final clean**********\n");
    }


    pre_clean(sources, coverage_cut, sg, bubble_dist);


    asg_arc_del_short_diploi_by_suspect_edge(sg, asm_opt.max_short_tip);
    asg_cut_tip(sg, asm_opt.max_short_tip);
    asg_arc_del_triangular_directly(sg, asm_opt.max_short_tip, reverse_sources, ruIndex);


    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_orthology_multiple_way(sg, reverse_sources, 0.4, asm_opt.max_short_tip, ruIndex);
    asg_cut_tip(sg, asm_opt.max_short_tip);

    



    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_too_short_overlaps(sg, 2000, min_ovlp_drop_ratio, reverse_sources, 
    asm_opt.max_short_tip, ruIndex);
    asg_cut_tip(sg, asm_opt.max_short_tip);

    asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, 100, 0);

    if (asm_opt.flag & HA_F_VERBOSE_GFA)
    {
        /*******************************for debug***************************************/
        write_debug_graph(sg, sources, coverage_cut, output_file_name, n_read, reverse_sources, ruIndex);
        debug_gfa:;
        /*******************************for debug***************************************/
    }
    
    /**
    debug_ma_hit_t(sources, coverage_cut, n_read, max_hang_length, 
    mini_overlap_length);
    debug_ma_hit_t(reverse_sources, coverage_cut, n_read, max_hang_length, 
    mini_overlap_length);
    **/
 
    ///note: don't apply asg_arc_del_too_short_overlaps() after this function!!!!
    rescue_contained_reads_aggressive(NULL, sg, sources, coverage_cut, ruIndex, max_hang_length, 
    mini_overlap_length, bubble_dist, 10, 1, 0, NULL, NULL);
    rescue_missing_overlaps_aggressive(NULL, sg, sources, coverage_cut, ruIndex, max_hang_length,
    mini_overlap_length, bubble_dist, 1, 0, NULL);
    rescue_missing_overlaps_backward(NULL, sg, sources, coverage_cut, ruIndex, max_hang_length, 
    mini_overlap_length, bubble_dist, 10, 1, 0);
    // rescue_wrong_overlaps_to_unitigs(NULL, sg, sources, reverse_sources, coverage_cut, ruIndex, 
    // max_hang_length, mini_overlap_length, bubble_dist, NULL);
    // rescue_no_coverage_aggressive(sg, sources, reverse_sources, &coverage_cut, ruIndex, max_hang_length, 
    // mini_overlap_length, bubble_dist, 10);

    if (ha_opt_triobin(&asm_opt))
    {
		char *buf = (char*)calloc(strlen(output_file_name) + 25, 1);
		sprintf(buf, "%s.dip", output_file_name);
        output_unitig_graph(sg, coverage_cut, buf, sources, ruIndex, max_hang_length, mini_overlap_length);
		free(buf);
        
        output_trio_unitig_graph(sg, coverage_cut, output_file_name, FATHER, sources, 
        reverse_sources, bubble_dist, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex, 
        0.05, 0.9, max_hang_length, mini_overlap_length);
        output_trio_unitig_graph(sg, coverage_cut, output_file_name, MOTHER, sources,
        reverse_sources, bubble_dist, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex, 
        0.05, 0.9, max_hang_length, mini_overlap_length);
    }
    else
    {

        output_unitig_graph(sg, coverage_cut, output_file_name, sources, ruIndex, max_hang_length, mini_overlap_length);

        if(VERBOSE >= 1)
        {
            output_read_graph(sg, coverage_cut, output_file_name, n_read);
        }

        output_contig_graph_primary_pre(sg, coverage_cut, output_file_name, sources, reverse_sources, 
        asm_opt.small_pop_bubble_size, asm_opt.max_short_tip, ruIndex, max_hang_length, mini_overlap_length);

        output_contig_graph_primary(sg, coverage_cut, output_file_name, sources, reverse_sources, 
        bubble_dist, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex, 0.05, 0.9, max_hang_length,
        mini_overlap_length);

        output_contig_graph_alternative(sg, coverage_cut, output_file_name, sources, ruIndex, max_hang_length, mini_overlap_length);
    }

	*coverage_cut_ptr = coverage_cut;
	*sg_ptr = sg;
}


void hamt_clean_graph(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, long long mini_overlap_length, 
long long max_hang_length, long long clean_round, long long gap_fuzz,
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, char* output_file_name, 
long long bubble_dist, int read_graph, R_to_U* ruIndex, asg_t **sg_ptr,   // the **sg_ptr is a NULL pointer
ma_sub_t **coverage_cut_ptr, int debug_g)
{
	ma_sub_t *coverage_cut = *coverage_cut_ptr;
	asg_t *sg = *sg_ptr;
    ma_ug_t *hamt_ug = 0;
    int simplereport[5];
    int acc;
    int cleanID = 0;
    vu32_t *long_tigs_in_resuce;
    double T0;

    paf_ct_v *pafidx_cis, *pafidx_trans;
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

    // extra rescues and others
    if (asm_opt.is_mode_low_cov){
        hamt_smash_haplotype(sources, reverse_sources, n_read);
    }

    if (debug_g && (asm_opt.do_probe_gfa==1)) {
        fprintf(stderr, "[M::%s] jump to proble dbg gfa\n", __func__);
        goto probe;
    }
    if(debug_g && !asm_opt.write_new_graph_bins) goto debug_gfa;
    fprintf(stderr, "[M::%s] no debug gfa\n", __func__);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

    renew_graph_init(sources, reverse_sources, sg, coverage_cut, ruIndex, n_read);
    reset_paf_buffers_sorting_marks_targetID(sources, n_read);
    reset_paf_buffers_sorting_marks_targetID(reverse_sources, n_read);

    // generate a indexing for pafs, in order to avoid get_specific_overlap calls
    T0 = Get_T();
    pafidx_cis = hamt_index_pafs_multithread(sources, 0, coverage_cut, n_read, 
                                            2048, asm_opt.thread_num) ;
    pafidx_trans = hamt_index_pafs_multithread(0, reverse_sources, coverage_cut, n_read, 
                                            2048, asm_opt.thread_num) ;
    fprintf(stderr, "[M::%s] generated paf index, used %.2f s\n", __func__, Get_T()-T0);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

    T0 = Get_T();
    hamt_symmetrize_paf(sources, pafidx_cis, n_read, asm_opt.thread_num);
    hamt_symmetrize_paf(reverse_sources, pafidx_trans, n_read, asm_opt.thread_num);
    fprintf(stderr, "[M::%s] symmetrize used %.2f s total.\n", __func__, Get_T()-T0);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

    asm_opt.get_specific_overlap_is_use_bf = 0; 
    //sort_paf_buffers_by_targetID(sources, n_read); sort_paf_buffers_by_targetID(reverse_sources, n_read);
    //hamt_normalize_ma_hit_t_single_side_advance(sources, n_read);
    //hamt_normalize_ma_hit_t_single_side_advance(reverse_sources, n_read);
    hamt_normalize_paf(sources, pafidx_cis, n_read, asm_opt.thread_num);
    hamt_normalize_paf(reverse_sources, pafidx_trans, n_read, asm_opt.thread_num);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

    memset(R_INF.trio_flag, AMBIGU, R_INF.total_reads*sizeof(uint8_t));





    hamt_clean_weak_ma_hit_t2(sources, reverse_sources, coverage_cut, n_read, 
            pafidx_cis, pafidx_trans);  
    ma_hit_sub(min_dp, sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());   





    detect_chimeric_reads_conservative(sources, n_read, readLen, coverage_cut, asm_opt.max_ov_diff_final * 2.0);  // detect_chimeric_reads(sources, n_read, readLen, coverage_cut, asm_opt.max_ov_diff_final * 2.0);
    ma_hit_cut(sources, n_read, readLen, mini_overlap_length, &coverage_cut);
    ma_hit_flt(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);

    if (!asm_opt.no_containedreads_heuristics){
        asm_opt.get_specific_overlap_is_use_bf = 0;
        hamt_hit_contained_multi(sources, reverse_sources, n_read, readLen, coverage_cut);
        hamt_hit_contained_drop_singleton_multi(sources, reverse_sources, n_read, readLen, coverage_cut);
    } else {
        fprintf(stderr, "[M::%s] skipped conatined reads sparing heuristics.\n", __func__);
    }

    hamt_ma_hit_contained_advance(sources, n_read, coverage_cut, ruIndex, max_hang_length, mini_overlap_length);  
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());


    // cleanup unuseful index
    destroy_paf_ct_v(pafidx_cis);
    destroy_paf_ct_v(pafidx_trans);

    sg = ma_sg_gen(sources, n_read, coverage_cut, max_hang_length, mini_overlap_length);

    asg_arc_del_trans(sg, gap_fuzz);


    asm_opt.coverage = get_coverage(sources, coverage_cut, n_read);

    if(VERBOSE >= 1)
    {
        char* unlean_name = (char*)malloc(strlen(output_file_name)+25);
        sprintf(unlean_name, "%s.unclean", output_file_name);
        output_read_graph(sg, coverage_cut, unlean_name, n_read);
        output_unitig_graph(sg, coverage_cut, unlean_name, sources, ruIndex, max_hang_length, mini_overlap_length);
        free(unlean_name);
    }

    asg_cut_tip(sg, asm_opt.max_short_tip);
    
    fprintf(stderr, "[M::%s] ====== initial clean ======\n", __func__);
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
    if(clean_round > 0)
    {
        double cut_step;
        if(clean_round == 1)
        {
            cut_step = max_ovlp_drop_ratio;
        }
        else
        {
            cut_step = (max_ovlp_drop_ratio - min_ovlp_drop_ratio) / (clean_round - 1);
        }
        double drop_ratio = min_ovlp_drop_ratio;
        int i = 0;
        for (i = 0; i < clean_round; i++, drop_ratio += cut_step)  // main clean loop
        {
            if(drop_ratio > max_ovlp_drop_ratio)
            {
                drop_ratio = max_ovlp_drop_ratio;
            }

            fprintf(stderr, "\n\n**********%d-th round drop: drop_ratio = %f**********\n", 
                            i, drop_ratio);

            hamt_asgarc_drop_tips_and_bubbles(sources, sg, 3, -1);  // pre_clean(sources, coverage_cut, sg, bubble_dist);
            
            /****************************may have bugs********************************/
            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_false_node_meta(sg, asm_opt.max_short_tip);
            asg_cut_tip(sg, asm_opt.max_short_tip);

            asg_arc_identify_simple_bubbles_multi(sg, 0);
            asg_arc_del_short_diploid_by_exact(sg, asm_opt.max_short_tip, sources);
            asg_cut_tip(sg, asm_opt.max_short_tip);
            /****************************may have bugs********************************/

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_diploid_by_length(sg, drop_ratio, asm_opt.max_short_tip, reverse_sources, 
            asm_opt.max_short_tip, 1, 1, 0, 0, ruIndex);

            asg_cut_tip(sg, asm_opt.max_short_tip);

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_short_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, 
            asm_opt.max_short_tip, ruIndex);     

            asg_arc_identify_simple_bubbles_multi(sg, 1);
            asg_arc_del_complex_false_link(sg, 0.6, 0.85, bubble_dist, reverse_sources, asm_opt.max_short_tip);

            asg_cut_tip(sg, asm_opt.max_short_tip);

        }
    }

    fprintf(stderr, "\n\n********** last round **********\n");

    hamt_asgarc_drop_tips_and_bubbles(sources, sg, 3, -1);  // pre_clean(sources, coverage_cut, sg, bubble_dist);

    asg_arc_del_short_diploi_by_suspect_edge(sg, asm_opt.max_short_tip);
    asg_cut_tip(sg, asm_opt.max_short_tip);
    asg_arc_del_triangular_directly(sg, asm_opt.max_short_tip, reverse_sources, ruIndex);
    asg_arc_identify_simple_bubbles_multi(sg, 0);
    asg_arc_del_orthology_multiple_way(sg, reverse_sources, 0.4, asm_opt.max_short_tip, ruIndex);

    asg_cut_tip(sg, asm_opt.max_short_tip);
    asg_arc_identify_simple_bubbles_multi(sg, 0);

    asg_arc_del_too_short_overlaps(sg, 2000, min_ovlp_drop_ratio, reverse_sources, 
                                   asm_opt.max_short_tip, ruIndex);

    asg_cut_tip(sg, asm_opt.max_short_tip);

    asg_arc_del_simple_circle_untig(sources, coverage_cut, sg, 100, 0);

    fprintf(stderr, "\n\n********** checkpoint: r_utg **********\n");
    fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
    if ((asm_opt.flag & HA_F_VERBOSE_GFA) || asm_opt.write_new_graph_bins)
    {
        /*******************************for debug***************************************/
        write_debug_graph(sg, sources, coverage_cut, output_file_name, n_read, reverse_sources, ruIndex);
        debug_gfa:;
        /*******************************for debug***************************************/
    }
    
    /**
    debug_ma_hit_t(sources, coverage_cut, n_read, max_hang_length, 
    mini_overlap_length);
    debug_ma_hit_t(reverse_sources, coverage_cut, n_read, max_hang_length, 
    mini_overlap_length);
    **/
 
    ///note: don't apply asg_arc_del_too_short_overlaps() after this function!!!!


    {
        // r_utg
        output_unitig_graph(sg, coverage_cut, output_file_name, sources, ruIndex, max_hang_length, mini_overlap_length);
        if(VERBOSE >= 1){
            output_read_graph(sg, coverage_cut, output_file_name, n_read);
        }

        // clean up ug
        {
            double time;
            hamt_asg_reset_seq_label(sg, 0);


            hamt_ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);  // lazy way of initing the trailing subgraph IDs

            // hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "HEY", cleanID);
            // hamt_ug_resolveTangles_threaded(sg, hamt_ug, asm_opt.thread_num, 0);
            // hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            // hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "HEYafter", cleanID);
            // exit(0);

            // topo pre clean
            fprintf(stderr, "[M::%s] ======= preclean =======\n", __func__);
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            acc = 0;
            for (int i=0; i<10; i++){
                time = Get_T();
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "before_initTopo_clnA", cleanID);}
                
                acc = hamt_ug_basic_topoclean(sg, hamt_ug, 0, 1, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "before_initTopo_clnB", cleanID);}

                // {  // debug
                //     kvec_asg_arc_t_warp new_rtg_edges;
                //     kv_init(new_rtg_edges.a);
                //     ma_ug_seq(hamt_ug, sg, &R_INF, coverage_cut, 
                //       sources, &new_rtg_edges, max_hang_length, mini_overlap_length);
                //     hamt_ug_opportunistic_elementary_circuits(sg, hamt_ug);
                //     exit(0);
                // }

                acc += hamt_ug_drop_redundant_nodes(sg, hamt_ug, 100000, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "before_initTopo_clnC", cleanID); cleanID++;}

                if (asm_opt.is_aggressive){
                    acc += hamt_ug_drop_transitive(sg, hamt_ug, 100000, 0);
                    hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                }

                acc += hamt_ug_cut_very_short_multi_tip(sg, hamt_ug, 4);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] round %d, dropped %d, used %.1fs\n", __func__, i, acc, Get_T()-time);
                if (acc==0) break;
            }
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_initTopo_cln", cleanID); cleanID++;}

            
            fprintf(stderr, "\n\n********** checkpoint: p_utg **********\n\n");
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            hamt_output_unitig_graph_advance(sg, coverage_cut, asm_opt.output_file_name, "p_utg", "utg",
                                     sources, ruIndex, max_hang_length, mini_overlap_length, -1);

            time = Get_T();
            // cut dangling circles and inversion links
            hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "HEYbefore_circle_cln", cleanID);
            hamt_circle_cleaning(sg, hamt_ug, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_circle_cln", cleanID); cleanID++;}

            // hap resuce
            hamt_ug_rescue_bifurTip(sg, hamt_ug, 0, sources, reverse_sources, coverage_cut);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_rescueBifurTip", cleanID); cleanID++;}

            hamt_ug_drop_shorter_ovlp(sg, hamt_ug, sources, reverse_sources);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);

            hamt_asgarc_ugCovCutDFSCircle_aggressive(sg, hamt_ug, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_covcutDFS_cln", cleanID); cleanID++;}

            fprintf(stderr, "\n\n[M::%s] time check #1, %.1fs\n", __func__, Get_T()-time);

            // topo
            time = Get_T();
            hamt_ug_basic_topoclean_simple(sg, hamt_ug, 0, 1, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_drop_redundant_nodes(sg, hamt_ug, 100000, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_TOPO2", cleanID); cleanID++;}

            hamt_ug_drop_redundant_nodes_bruteforce(sg, hamt_ug, 100000, 0, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_TOPO2-brute", cleanID); cleanID++;}
            fprintf(stderr, "\n\n[M::%s] time check #2, %.1fs\n", __func__, Get_T()-time);

            // more topo cleaning
            time = Get_T();
            hamt_ug_prectgTopoClean(sg, coverage_cut, sources, ruIndex, 0, 1, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "after_TOPO2_and_prectg", cleanID); cleanID++;}
            hamt_ug_oneutgCircleCut(sg, hamt_ug, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "TOPO3_before", cleanID); cleanID++;}
            hamt_ug_basic_topoclean_simple(sg, hamt_ug, 0, 1, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_drop_redundant_nodes(sg, hamt_ug, 100000, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_drop_redundant_nodes_bruteforce(sg, hamt_ug, 100000, 0, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "\n\n[M::%s] time check #3, %.1fs\n", __func__, Get_T()-time);

            // resolve complex bubble
            fprintf(stderr, "\n\n[M::%s] ======= complex bubbles =======\n", __func__);
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            int nb_complex_bubble_cut;
            for (int round_resolve=0; round_resolve<5; round_resolve++){
                time = Get_T();
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "resolveCmplx_before", cleanID); cleanID++;}
                
                nb_complex_bubble_cut = hamt_ug_prectg_resolve_complex_bubble(sg, hamt_ug, 0, 1, 0, asm_opt.gc_superbubble_tig_max_length);
                nb_complex_bubble_cut += hamt_ug_drop_midsizeTips(sg, hamt_ug, 5, 0);
                nb_complex_bubble_cut += hamt_ug_drop_midsizeTips_aggressive(sg, hamt_ug, 0.7, 0);
                nb_complex_bubble_cut += hamt_ug_resolve_oneMultiLeafSoapBubble(sg, hamt_ug, 0, 1, 0);
                nb_complex_bubble_cut += hamt_ug_resolve_small_multileaf_with_covcut(sg, hamt_ug, 75000, 2, 0);
                nb_complex_bubble_cut += hamt_ug_pop_simpleInvertBubble(sg, hamt_ug, 0, 1, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                
                hamt_ug_drop_redundant_nodes(sg, hamt_ug, 200000, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);

                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "resolveCmplx_before-brute", cleanID); cleanID++;}
                hamt_ug_drop_redundant_nodes_bruteforce(sg, hamt_ug, 200000, 0, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
               
                fprintf(stderr, "[M::%s] dropped %d, used %.1fs\n", __func__, nb_complex_bubble_cut, Get_T()-time);

                if (nb_complex_bubble_cut==0){
                    break;
                }
            }

            fprintf(stderr, "\n\n[M::%s] ======= circle cut and topo cleans =======\n", __func__);
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            time = Get_T();
            {
                int tot=0;
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                tot += hamt_ug_oneutgCircleCut2(sg, hamt_ug, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                tot += hamt_ug_oneutgCircleCut(sg, hamt_ug, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                tot += hamt_ug_special_toploclean(sg, hamt_ug, 0, 1);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] treated %d, used %.1fs\n", __func__, tot, Get_T()-time);
            }
            if (asm_opt.do_probe_gfa<0){  // --probe-gfa specified but bin files do no yet exist
                hamt_write_debug_graph(sg, sources, coverage_cut, output_file_name, n_read, reverse_sources, ruIndex);
            }
probe:
            if (!hamt_ug) hamt_ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            // resolve tangle
            fprintf(stderr, "\n\n[M::%s] ======= tangles =======\n", __func__);
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            int nb_tangle_cut, nb_tip;
            for (int round_resolve=0; round_resolve<5; ){
                time = Get_T();
                round_resolve++;
                /*if (asm_opt.write_debug_gfa)*/ {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "resolveTangle_before", cleanID); cleanID++;}
                
                time = Get_T();
                nb_tangle_cut = hamt_ug_resolveTangles_threaded(sg, hamt_ug, asm_opt.thread_num, 0, round_resolve);
                if (nb_tangle_cut>0) {
                    round_resolve--;
                    fprintf(stderr, "[M::%s] segment #1: will extend one round\n", __func__);
                }
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] segment #1: %.1fs\n", __func__, Get_T()-time);

                time = Get_T();
                nb_tangle_cut = hamt_ug_cut_shortTips_arbitrary(sg, hamt_ug, 50000, 0);
                fprintf(stderr, "[M::%s] segment #2: %.1fs\n", __func__, Get_T()-time);
                
                time = Get_T();
                nb_tangle_cut += hamt_ug_covcut_falseVertexLoop(sg, hamt_ug, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] segment #3: %.1fs\n", __func__, Get_T()-time);

                time = Get_T();
                hamt_ug_drop_redundant_nodes(sg, hamt_ug, 200000, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] segment #4: %.1fs\n", __func__, Get_T()-time);

                time = Get_T();
                hamt_ug_drop_redundant_nodes_bruteforce(sg, hamt_ug, 200000, 0, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                hamt_circle_cleaning(sg, hamt_ug, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] segment #5: %.1fs\n", __func__, Get_T()-time);

                time = Get_T();
                nb_tangle_cut += hamt_ug_disconnect_long_contig_pairs_by_cov(sg, hamt_ug); 
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                fprintf(stderr, "[M::%s] segment #6: %.1fs\n", __func__, Get_T()-time);

                fprintf(stderr, "[M::%s] round %d, cut %d, used %.1fs\n", __func__, round_resolve, nb_tangle_cut, Get_T()-time);
                if (nb_tangle_cut==0){
                    break;
                }
            }

            // rescues
            fprintf(stderr, "\n\n[M::%s] ======= rescues and coverage-based cut =======\n", __func__);
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            time = Get_T();
            {
                int tot=0;
                tot += hamt_ug_pop_simpleShortCut(sg, hamt_ug, 0, 1, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                
                tot += hamt_ug_prectg_rescueShortCircuit(sg, sources, reverse_sources, ruIndex, coverage_cut, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "afterRescue", cleanID);cleanID++;}
                
                tot += hamt_ug_prectg_rescueShortCircuit_simpleAggressive(sg, hamt_ug, sources, reverse_sources, coverage_cut, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
                if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "afterRescue", cleanID);cleanID++;}
                
                fprintf(stderr, "[M::%s] treated %d, used %.1fs\n", __func__, tot, Get_T()-time);
            }

            // hap cov cut
            time = Get_T();
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "beforeHapCovCut", cleanID);cleanID++;}
            hamt_ug_treatBifurcation_hapCovCut(sg, hamt_ug, 0.7, 0.5, reverse_sources, 0, 1);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "midHapCovCut", cleanID);cleanID++;}
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "[M::%s] time check #1: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_try_circularize(sg, hamt_ug, sources, reverse_sources, ruIndex, coverage_cut, 1000000);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "[M::%s] time check #2: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ugasg_cut_shortTips(sg, hamt_ug, 0, 1, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_cut_shortTips_arbitrary(sg, hamt_ug, 30000, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "[M::%s] time check #3: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "beforeBIF", cleanID);cleanID++;}
            hamt_ug_resolve_fake_haplotype_bifurcation(sg, hamt_ug, 0, sources, reverse_sources);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "beforeBIFagg", cleanID); cleanID++;}
            fprintf(stderr, "[M::%s] time check #4: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_pop_tinyFlatCircles(sg, hamt_ug, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "afterflatcircle", cleanID); cleanID++;}
            fprintf(stderr, "[M::%s] time check #5: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_prectg_rescueShortCircuit(sg, sources, reverse_sources, ruIndex, coverage_cut, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "afterShortCircuit", cleanID); cleanID++;}
            fprintf(stderr, "[M::%s] time check #6: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_pop_unevenInvertBubble(sg, hamt_ug, 0, 1);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(hamt_ug, asm_opt.output_file_name, "afterunevenInv", cleanID); cleanID++;}
            fprintf(stderr, "[M::%s] time check #7: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_cleanup_almost_circular(sg, hamt_ug, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_pop_simpleInvertBubble(sg, hamt_ug, 0, 1, 0);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "[M::%s] time check #8: used %.1fs\n", __func__, Get_T()-time);

            time = Get_T();
            hamt_ug_try_circularize(sg, hamt_ug, sources, reverse_sources, ruIndex, coverage_cut, 1000000);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            hamt_ug_rescueLongUtg(sg, sources, reverse_sources, ruIndex, coverage_cut);
            hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
            fprintf(stderr, "[M::%s] time check #9: used %.1fs\n", __func__, Get_T()-time);
        }  

        fprintf(stderr, "\n\n********** checkpoint: p_ctg **********\n\n");
        fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());

        // p_ctg and a_ctg
        if (!asm_opt.is_use_exp_graph_cleaning){
            output_contig_graph_primary(sg, coverage_cut, output_file_name, sources, reverse_sources, 
            bubble_dist, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex, 0.05, 0.9, max_hang_length,
            mini_overlap_length);

            output_contig_graph_alternative(sg, coverage_cut, output_file_name, sources, ruIndex, max_hang_length, mini_overlap_length);
        }else{

            // contig graph
            simplereport[0] = hamt_ug_delete_unconnected_single_read_contigs(sg, hamt_ug);
            hamt_ug_destroy(hamt_ug);
            hamt_ug = hamt_ug_gen(sg, coverage_cut, reverse_sources, ruIndex, 0);
            simplereport[1] = hamt_ug_delete_unconnected_single_read_contigs(sg, hamt_ug);
            fprintf(stderr, "[M::%s] removed single-read unconnected contigs, %d on primary, %d on alt.\n", __func__, simplereport[0], simplereport[1]);
            
            hamt_output_unitig_graph_advance(sg, coverage_cut, asm_opt.output_file_name, "p_ctg", "ctg",
                                     sources, ruIndex, max_hang_length, mini_overlap_length, 0);
            hamt_output_unitig_graph_advance(sg, coverage_cut, asm_opt.output_file_name, "a_ctg", "ctg",
                                     sources, ruIndex, max_hang_length, mini_overlap_length, 1);
            //if (asm_opt.do_probe_gfa<0){  // --probe-gfa specified but bin files do no yet exist
            //    hamt_write_debug_graph(sg, sources, coverage_cut, output_file_name, n_read, reverse_sources, ruIndex);
            //}
            
            fprintf(stderr, "\n\n********** checkpoint: post-assembly **********\n\n");
//probe:
            fprintf(stderr, "[M::%s] (peak RSS so far: %.1f GB)\n", __func__, Get_U());
            if (!asm_opt.no_post_assembly_binning){
                if (!hamt_ug) hamt_ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, 0);
                hamt_ug_regen(sg, &hamt_ug, coverage_cut, sources, ruIndex, 0);
    
                // path finding-based circle rescue
                kvec_asg_arc_t_warp new_rtg_edges;
                kv_init(new_rtg_edges.a);
                ma_ug_seq(hamt_ug, sg, &R_INF, coverage_cut, 
                          sources, &new_rtg_edges, max_hang_length, mini_overlap_length);
                long_tigs_in_resuce = hamt_ug_opportunistic_elementary_circuits(sg, hamt_ug, asm_opt.thread_num);
                kv_destroy(new_rtg_edges.a);
                
                // simple binning
                hamt_simple_binning(hamt_ug, long_tigs_in_resuce, asm_opt.thread_num, 
                        asm_opt.output_file_name, asm_opt.write_binning_fasta);
                kv_destroy(*long_tigs_in_resuce);
            }else{
                fprintf(stderr, "[M::%s] skipped binning. \n", __func__);
            }


        if (hamt_ug) {hamt_ug_destroy(hamt_ug);}
        else {fprintf(stderr, "WARNING: ug did not exist, were you using --probe-gfa?\n");}
        }
    }


	*coverage_cut_ptr = coverage_cut;
	*sg_ptr = sg;

    if (!asm_opt.is_use_exp_graph_cleaning) {fprintf(stderr, "Inconsistency threshold for low-quality regions in BED files: %u%%\n", asm_opt.bed_inconsist_rate);}
}




void build_string_graph_without_clean(
int min_dp, ma_hit_t_alloc* sources, ma_hit_t_alloc* reverse_sources, 
long long n_read, uint64_t* readLen, long long mini_overlap_length, 
long long max_hang_length, long long clean_round, long long gap_fuzz,
float min_ovlp_drop_ratio, float max_ovlp_drop_ratio, char* output_file_name, 
long long bubble_dist, int read_graph, int write)
{
    R_to_U ruIndex;
    init_R_to_U(&ruIndex, n_read);
    asg_t *sg = NULL;
    ma_sub_t* coverage_cut = NULL;
    int loaded_debug_gfa;
        
    ///actually min_thres = asm_opt.max_short_tip + 1 there are asm_opt.max_short_tip reads
    min_thres = asm_opt.max_short_tip + 1;

    if ( (asm_opt.flag & HA_F_VERBOSE_GFA) || asm_opt.do_probe_gfa)
    {
        fprintf(stderr, "try loading bin files of assembly graph...\n");
        fflush(stderr);
        if(asm_opt.do_probe_gfa){
            loaded_debug_gfa = hamt_load_debug_graph(&sg, &sources, &coverage_cut, asm_opt.probebin_base_name, &reverse_sources, &ruIndex, n_read);    
            if (!loaded_debug_gfa){  // fallback
                fprintf(stderr, "--probe-gfa loading failed, fallback to try --dbg-gfa...\n");
                asm_opt.do_probe_gfa = -1;
                loaded_debug_gfa = load_debug_graph(&sg, &sources, &coverage_cut, asm_opt.bin_base_name, &reverse_sources, &ruIndex, n_read);
            }
        }else{
            double T = Get_T();
            loaded_debug_gfa = load_debug_graph(&sg, &sources, &coverage_cut, asm_opt.bin_base_name, &reverse_sources, &ruIndex, n_read);
            fprintf(stderr, "[M::%s] load_debug_graph(--dbg-gfa) used %.2f s\n", 
                    __func__, Get_T()-T); fflush(stderr);
        }
        if(loaded_debug_gfa)
        {
            fprintf(stderr, "debug gfa has been loaded (code %d)\n", loaded_debug_gfa);
            
            if (!asm_opt.is_use_exp_graph_cleaning){
                clean_graph(min_dp, sources, reverse_sources, n_read, readLen, mini_overlap_length, 
                max_hang_length, clean_round, gap_fuzz, min_ovlp_drop_ratio, max_ovlp_drop_ratio, 
                output_file_name, bubble_dist, read_graph, &ruIndex, &sg, &coverage_cut, 1);
            }else{
                hamt_clean_graph(min_dp, sources, reverse_sources, n_read, readLen, mini_overlap_length, 
                max_hang_length, clean_round, gap_fuzz, min_ovlp_drop_ratio, max_ovlp_drop_ratio, 
                output_file_name, bubble_dist, read_graph, &ruIndex, &sg, &coverage_cut, 1);
            }
            asg_destroy(sg);
            free(coverage_cut);
            destory_R_to_U(&ruIndex);
            return;
        }
    }
    
    // output stuff before generate graph
    if (asm_opt.write_index_to_disk && write)  // write_index_to_disk is always 1; write is !overlap_loaded
    {
        write_all_data_to_disk(sources, reverse_sources, 
        &R_INF, output_file_name);
    }

    ///debug_info_of_specfic_read("m64062_190803_042216/177341795/ccs", sources, reverse_sources, -1, "beg");
    // debug_info_of_specfic_read("m64062_190807_194840/126682874/ccs", sources, reverse_sources, -1, "beg");

    if (!(asm_opt.flag & HA_F_BAN_ASSEMBLY))
    {
        sort_paf_buffers_by_targetID(sources, n_read);
        sort_paf_buffers_by_targetID(reverse_sources, n_read);
        try_rescue_overlaps(sources, reverse_sources, n_read, 4); 

        if (!asm_opt.is_use_exp_graph_cleaning){
            clean_graph(min_dp, sources, reverse_sources, n_read, readLen, mini_overlap_length, 
            max_hang_length, clean_round, gap_fuzz, min_ovlp_drop_ratio, max_ovlp_drop_ratio, 
            output_file_name, bubble_dist, read_graph, &ruIndex, &sg, &coverage_cut, 0);
        }else{
            hamt_clean_graph(min_dp, sources, reverse_sources, n_read, readLen, mini_overlap_length, 
            max_hang_length, clean_round, gap_fuzz, min_ovlp_drop_ratio, max_ovlp_drop_ratio, 
            output_file_name, bubble_dist, read_graph, &ruIndex, &sg, &coverage_cut, 0);
        }


        asg_destroy(sg);
        free(coverage_cut);
    }

    destory_R_to_U(&ruIndex);
}

