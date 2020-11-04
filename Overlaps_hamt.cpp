#include <stdint.h>
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>
#include <assert.h>
#include <pthread.h>
#include "Overlaps.h"
#include "Process_Read.h"
#include "CommandLines.h" 
#include "Overlaps_hamt.h"
#include "ksort.h"

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)
KRADIX_SORT_INIT(ovhamt64, uint64_t, uint64_t, 8)
KRADIX_SORT_INIT(ovhamt32, uint32_t, uint32_t, 4)

//////////////////////////////////////////////////////////////////////
//                        debug functions                           //
//////////////////////////////////////////////////////////////////////
void write_debug_auxsg(asg_t *sg, char *file_base_name){
    char* index_name = (char*)malloc(strlen(file_base_name)+15);
    sprintf(index_name, "%s.auxsg.gfa", file_base_name);
    FILE* fp = fopen(index_name, "w");

    uint32_t v, w, dir1, dir2;
    for (uint32_t i=0; i<sg->n_arc; i++){
        v = sg->arc[i].ul>>32;
        w = sg->arc[i].v;
        dir1 = (uint32_t)v&1;
        dir2 = (uint32_t)w&2;
        fprintf(fp, "#%d\t%d\t#%d\t%d\n", (int)v+1, (int)dir1, (int)w+1, (int)dir2);
    }
    fclose(fp);
    free(index_name);
}

void write_debug_assembly_graph(asg_t *sg, All_reads *rs, char* read_file_name){
    fprintf(stderr, "Writing debug asg to disk... \n");
    double startTime = Get_T();
    char* index_name = (char*)malloc(strlen(read_file_name)+15);
    sprintf(index_name, "%s.dbg_asg", read_file_name);
    FILE* fp = fopen(index_name, "w");

    // asg64_v a = {0,0,0};
	uint32_t n_vtx = sg->n_seq * 2, v, i;
    
    uint32_t nv;
    asg_arc_t* av;

    // basic info of seqs
    for (i=0; i<rs->total_reads; i++){
        if (sg->seq[i].del) continue;
        fprintf(fp, "s\t%" PRIu32 "\t%f\t%" PRIu16 "\t%f\t%" PRIu8 "\t%" PRIu8 "\n",
                           i, rs->mean[i], rs->median[i], rs->std[i], rs->mask_readtype[i], rs->mask_readnorm[i]);
    }
    
    // the graph
	for (v = 0; v < n_vtx; ++v) {
        if (sg->seq[v>>1].del) continue;
        nv = asg_arc_n(sg, v);
	    av = asg_arc_a(sg, v);
        if (nv==0) continue;
        for (i=0; i<nv; i++){
            if (av[i].del) continue;
            if (sg->seq[av[i].v>>1].del) continue;
            fprintf(fp, "l\t%" PRIu32 "\t%d\t%" PRIu32 "\t%d\t", v>>1, v&1, av[i].v>>1, av[i].v&1);  // vid, v_strand, wid, w_strand
            fprintf(fp, "%" PRIu8 "\n", sg->seq_vis[i]);
        }
    }

    fflush(fp);
    fclose(fp);
    fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    free(index_name);
}

void write_debug_SCC_labeling(asg_t *g, int *SCC_labels, char *filename){
    // usg
    fprintf(stderr, "Writing SCC labels to disk...\n");
    char *filename_out = (char*)malloc(strlen(filename)+25);
    sprintf(filename_out, "%s.SCClabels", filename);
    
    FILE *fp = fopen(filename_out, "w");
    int SCC_label;
    uint32_t idx;
    // for (int i=0; i<g->n_arc; i++){
    for (int i=0; i<g->n_seq*2; i++){
        if (g->whitelist){
            if (!g->whitelist[i]){
                continue;
            }
        }
        idx = (uint32_t)i>>1;
        if (SCC_labels[idx]==0){
            continue;
            // fprintf(fp, "%.*s\t%c\torphan\n", (int)Get_NAME_LENGTH((R_INF), idx), Get_NAME((R_INF), idx),
            //                         "+-"[g->arc[i].ul>>32&1]);    
        }else{
            fprintf(fp, "%.*s\t%d\n", (int)Get_NAME_LENGTH((R_INF), idx), Get_NAME((R_INF), idx), SCC_labels[idx]);
        }
    }
    fclose(fp);
    
    free(filename_out);
    fprintf(stderr, "Wrote SCC labels.\n");
}

void write_debug_SCC_labeling2(ma_ug_t *ug, asg_t *g, int *SCC_labels, char *filename){
    // asg
    fprintf(stderr, "Writing SCC labels to disk...\n");
    char *filename_out = (char*)malloc(strlen(filename)+50);
    sprintf(filename_out, "%s.SCClabels_usgWithUtgName", filename);
    
    uint32_t *v2vu = (uint32_t*)calloc(g->n_seq*2, sizeof(uint32_t));
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        v2vu[ug->u.a[vu].start>>1] = vu+1;
        v2vu[ug->u.a[vu].end>>1] = vu+1;
    }

    FILE *fp = fopen(filename_out, "w");
    int SCC_label;
    uint32_t idx;
    for (int i=0; i<g->n_seq*2; i++){
        if (g->whitelist){
            if (!g->whitelist[i]){
                continue;
            }
        }
        idx = (uint32_t)i>>1;
        if (SCC_labels[i]==0){
            fprintf(fp, "what?\n");
            continue;
            // fprintf(fp, "%.*s\t%c\torphan\n", (int)Get_NAME_LENGTH((R_INF), idx), Get_NAME((R_INF), idx),
            //                         "+-"[g->arc[i].ul>>32&1]);    
        }else{
            fprintf(fp, "utg%.6d\t%.*s\t%d\n", v2vu[idx],(int)Get_NAME_LENGTH((R_INF), idx), Get_NAME((R_INF), idx), SCC_labels[i]);
        }
    }
    fclose(fp);
    
    free(filename_out);
    fprintf(stderr, "Wrote SCC labels.\n");
    free(v2vu);
}

void write_debug_SCC_labeling_ug(ma_ug_t *ug, int *SCC_labels, char *filename){
    fprintf(stderr, "[M::%s]Writing SCC labels to disk...\n", __func__);
    char *filename_out = (char*)malloc(strlen(filename)+25);
    sprintf(filename_out, "%s.SCClabels_rutg", filename);
    
    FILE *fp = fopen(filename_out, "w");
    uint32_t idx;
    char name[32];
    fprintf(fp, "#utgID\tSCC_label\n");
    for (int i=0; i<(int)ug->u.n; i++){
        sprintf(name, "utg%.6d%c", i + 1, "lc"[ug->u.a[i].circ]);
        fprintf(fp, "%s\t%d\n", name, SCC_labels[i]);
    }

    fclose(fp);
    free(filename_out);
    fprintf(stderr, "[M::%s]Wrote SCC labels.\n", __func__);
}


void ma_ug_print2_lite(const ma_ug_t *ug, All_reads *RNF, asg_t* read_g,
                       int print_seq, const char* prefix, FILE *fp)
{
    // lite: no seq, no coverage
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
        sprintf(name, "%s%.6d%c", prefix, i + 1, "lc"[p->circ]);
		if (print_seq) fprintf(fp, "S\t%s\t%s\tLN:i:%d\tdp:f:1\n", name, p->s? p->s : "*", p->len);
		else fprintf(fp, "S\t%s\t*\tLN:i:%d\tdp:f:1\n", name, p->len);
        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if(x<RNF->total_reads)
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
                (int)Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), 
                0, 10000, x, "apmaaa"[RNF->trio_flag[x]]);  // was coverage_cut[x].s and coverage_cut[x].e
            }
            else
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
					"FAKE", 0, 10000, x, '*');  // was coverage_cut[x].s and coverage_cut[x].e
            }
            
            
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
        prefix, (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
		prefix,	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}
}

void hamtdebug_add_fake_utg_coverage(ma_ug_t *ug){
    if (ug->utg_coverage){free(ug->utg_coverage);}
    ug->utg_coverage = (int*)calloc(ug->g->n_seq, sizeof(int));
}

void hamt_ug_print_simple(const ma_ug_t *ug, All_reads *RNF, FILE *fp)
{
    // debug function. Assume coverage has been collected. (otherwise it's ma_ug_print2)
	uint32_t i, j, l;
	char name[32];
	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
		ma_utg_t *p = &ug->u.a[i];
        if(p->m == 0) continue;
        sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
		fprintf(fp, "S\t%s\t*\tLN:i:%d\tdp:f:%u\n", name, p->len, ug->utg_coverage[i]);
        
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if(x<RNF->total_reads)
            {
                fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+-"[p->a[j]>>32&1],
                (int)Get_NAME_LENGTH((*RNF), x), Get_NAME((*RNF), x), 
                0, 8888, x, "apmaaa"[RNF->trio_flag[x]]);
            }else{
                fprintf(stderr, "[%s] ERROR\n", __func__);
                exit(1);
            }
        }
	}
	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\t%s%.6d%c\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
        "utg", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
		"utg",	(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}
}

void hamtdebug_output_unitig_graph_ug(ma_ug_t *ug, char *base_fn, int cleanID){
    // write gfa without ma_ug_gen
    
    // kvec_asg_arc_t_warp new_rtg_edges;
    // kv_init(new_rtg_edges.a);
    // ma_ug_seq(ug, sg, &R_INF, coverage_cut, sources, &new_rtg_edges, max_hang, min_ovlp);

    char* gfa_name = (char*)malloc(strlen(base_fn)+25);
    sprintf(gfa_name, "%s.G_%d.rutg.gfa", base_fn, cleanID);
    FILE *output_file = fopen(gfa_name, "w");
    hamt_ug_print_simple(ug, &R_INF, output_file);
    fclose(output_file);

    free(gfa_name);
    // kv_destroy(new_rtg_edges.a);
}

//////////////////////////////////////////////////////////////////////
//                        helper structs                           //
//////////////////////////////////////////////////////////////////////
typedef struct {
    int threadID, total_threads;
    asg_t *sg;
    ma_ug_t *ug;
} thrdpara_ugPopBubble;



//////////////////////////////////////////////////////////////////////
//                        helper routines                           //
//////////////////////////////////////////////////////////////////////


typedef struct {
    uint32_t *a;
    int n, m;  // n is the count, m is the capacity
}stacku32_t;

void stacku32_init(stacku32_t *stack){
    stack->n = 0;
    stack->m = 16;
    stack->a = (uint32_t*)malloc(16*sizeof(uint32_t));
    assert(stack->a);
}

void stacku32_destroy(stacku32_t *stack){
    free(stack->a);
}

void stacku32_push(stacku32_t *stack, uint32_t value){
    if (stack->n+2>stack->m){
        stack->m = stack->m<16? 16 : (stack->m + (stack->m>>1));
        stack->a = (uint32_t*)realloc(stack->a, sizeof(uint32_t)*stack->m);
        assert(stack->a);
    }
    stack->a[stack->n++] = value;
}

int stacku32_pop(stacku32_t *stack, uint32_t *value){
    if (stack->n==0){
        return 0;
    }
    *value = stack->a[stack->n-1];
    stack->n--;
    return 1;
}

void stacku32_reset(stacku32_t *stack){
    stack->n = 0;
}

typedef struct {
    uint32_t *a;
    int n, m;  // n is the count, m is the capacity
    int i_head, i_tail;
    int earlywrap, prev_m;  // because variable length
}queue32_t;

void queue32_init(queue32_t *q){
    q->a = (uint32_t*)malloc(sizeof(uint32_t)*16);
    assert(q->a);
    q->n = 0; 
    q->m = 16;
    q->i_head = 0;
    q->i_tail = 0;
    q->earlywrap = 0;
    q->prev_m = q->m;
}

void queue32_destroy(queue32_t *q){
    free(q->a);
}

void queue32_enqueue(queue32_t *q, uint32_t d){
    q->a[q->i_tail] = d;
    q->n++;
    if (q->i_tail==(q->m-1)){
        q->i_tail = 0;
    }else{
        q->i_tail++;
    }
    if ((q->n+2)>=q->m){
        q->a = (uint32_t*)realloc(q->a, sizeof(uint32_t)*(q->m + (q->m>>1)));
        assert(q->a);
        q->prev_m = q->m;
        q->m = q->m + (q->m>>1);
        if (q->i_head>q->i_tail){
            q->earlywrap = 1;
        }
    }
}

int queue32_dequeue(queue32_t *q, uint32_t *d){
    if (q->n==0){ // underflow
        return 0;
    }
    *d = q->a[q->i_head];
    q->n--;
    if (!q->earlywrap){  // regular case
        if (q->i_head==(q->m-1)){
            q->i_head = 0;
        }else{
            q->i_head++;
        }
    }else{
        if (q->i_head==(q->prev_m-1)){
            q->i_head = 0;
            q->earlywrap = 0;
        }else{
            q->i_head++;
        }
    }
    return 1;
}

void queue32_reset(queue32_t *q){
    q->n = 0;
    q->i_head = 0;
    q->i_tail = 0;
    q->earlywrap = 0;
    q->prev_m = q->m;
}

// set .del for the arc between unitig vu and wu (vu and wu both have direction, 
// which is as defined in ma_ug_gen a.k.a. as in ug->g)
static inline void ugasg_arc_del(asg_t *g, ma_ug_t *ug, uint32_t vu, uint32_t wu, int del)
{
    uint32_t v, w;
    if ((vu&1)==0){
        v = ug->u.a[vu>>1].end^1;
    }else{
        v = ug->u.a[vu>>1].start^1;
    }
    if ((wu&1)==0){
        w = ug->u.a[wu>>1].start;
    }else{
        w = ug->u.a[wu>>1].end;
    }
    asg_arc_del(g, v, w, del);
}

static inline void hamt_ug_arc_del(asg_t *sg, ma_ug_t *ug, uint32_t vu, uint32_t wu, int del){
    // drop arcs between vu and wu in both directions, on sg and ug
    asg_arc_del(ug->g, vu, wu, del);
    asg_arc_del(ug->g, wu^1, vu^1, del);
    ugasg_arc_del(sg, ug, vu, wu, del);
    ugasg_arc_del(sg, ug, wu^1, vu^1, del);
}

void hamt_collect_utg_coverage(asg_t *sg, ma_ug_t *ug, 
                                const ma_sub_t* coverage_cut,
                                ma_hit_t_alloc* sources, R_to_U* ruIndex){
    // FUNC
    //     collect unitig coverages and store in ug->utg_coverage
    double startTime = Get_T();

    if (ug->utg_coverage){
        free(ug->utg_coverage);
    }
    ug->utg_coverage = (int*)calloc(ug->u.n, sizeof(int));
    uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));
    for (uint32_t i=0; i<ug->u.n; i++){
        ug->utg_coverage[i] = get_ug_coverage(&ug->u.a[i], sg, coverage_cut, sources, ruIndex, primary_flag);
    }
    free(primary_flag);
    if (VERBOSE){
        fprintf(stderr, "[M::%s] collected ug coverages.\n", __func__);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}
void hamt_destroy_utg_coverage(ma_ug_t *ug){
    if (ug->utg_coverage){
        free(ug->utg_coverage);
        ug->utg_coverage = 0;
    }
}

void hamt_ug_regen(asg_t *sg, ma_ug_t **ug,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex){
    hamt_destroy_utg_coverage(*ug);
    ma_ug_destroy(*ug);
    *ug = ma_ug_gen(sg);
    hamt_collect_utg_coverage(sg, *ug, coverage_cut, sources, ruIndex);       
}



int hamt_check_covDiff(int cov_min, int cov_max){
    if (cov_max<=5){
        return 0;
    }
    if (cov_max<30){
        if ((cov_max-cov_min)>15 || ((float)cov_max/cov_min)>=2){
            return 1;
        }else{
            return 0;
        }
    }else if (cov_max<50){
        if ((cov_max-cov_min)>20 || ((float)cov_max/cov_min)>=1.5){
            return 1;
        }else{
            return 0;
        }
    }else{
        if ((cov_max-cov_min)>50 || ((float)cov_max/cov_min)>=1.5){
            return 1;
        }else{
            return 0;
        }
    }
}

uint32_t asg_arc_get_complementaryID(asg_t *sg, asg_arc_t *arc){
    uint32_t v = arc->v^1;
    uint32_t w = (arc->ul>>32)^1;
    uint32_t i;
    for (i=0; i<sg->n_arc; i++){
        if (((sg->arc[i].ul>>32)==v) && (sg->arc[i].v==w)){
            return i;
        }
    }
    fprintf(stderr, "[%s] unexpected.\n", __func__);
    exit(1);
}

/////////////////////////
//         asg         //
/////////////////////////

static inline int hamt_asgarc_util_countPre(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc){
    //?????? should i consider sg->seq[v].c as well?
    // count the predecessors of vertex v (i.e. the targets of v^1)
    if (include_del_seq && include_del_seq){
        return asg_arc_n(g, v^1);
    }
    if (!include_del_seq && g->seq[v>>1].del){
        return 0;
    }
    uint32_t nv, nv0 = asg_arc_n(g, v^1);
    asg_arc_t *av = asg_arc_a(g, v^1);
    int i;
    for (i = nv = 0; i<(int)nv0; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        nv++;
    }
    return nv;
}
static inline int hamt_asgarc_util_countSuc(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc){
    //?????? should i consider sg->seq[v].c as well?
    // count the targets of vertex v
    if (include_del_seq && include_del_seq){
        return asg_arc_n(g, v);
    }
    if (!include_del_seq && g->seq[v>>1].del){
        return 0;
    }
    uint32_t nv, nv0 = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    for (i = nv = 0; i<(int)nv0; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        nv++;
    }
    return nv;
}

// int hamt_asgarc_util_available_vtx_of_node(asg_t *g, uint32_t nodeID, uint32_t *v){
//     // given a node ID, store either of the vertices that has at least 1 arc to v
//     // return 1 if found, 0 otherwise
//     if (asg_arc_n(g, nodeID<<1)>0){
//         *v = nodeID<<1;
//         return 1;
//     }else if (asg_arc_n(g, nodeID<<1 | 1)>0){
//         *v = nodeID<<1 | 1;
//         return 1;
//     }else{
//         return 0;
//     }
// }


int hamt_asgutil_detect_strict_single_long_path(asg_t *sg, uint32_t begNode, uint32_t *endNode, 
                                                int threshold_nodes, uint64_t threshold_length){
    // simpler & more aggressive version of detect_single_path* functions.
    // - exclude any deleted seq or arc
    // - does not store any path info, only report its type
    // (set a threshold to zero will disable it)
    // returns: 0 is nope, 1 is yes, -1 is error
    if (threshold_nodes==0 && threshold_length==0 && ((!endNode) || ((*endNode)>>1 == begNode>>1))){
        fprintf(stderr, "[W::%s] no threshold/destination set.\n", __func__);
        return -1;
    }
    uint32_t v = begNode, nv;
    asg_arc_t *av;
    int cnt = 0;  // nodes visited
    uint64_t acc = 0;  // length of seqs visited
    uint32_t i;
    while (1){
        if (cnt>threshold_nodes || acc>threshold_length) break;
        if (endNode && (v>>1 == (*endNode)>>1)) break;
        nv = hamt_asgarc_util_countSuc(sg, v, 0, 0);
        if (nv>1) return 0;  // bifur
        if (nv==0) return 0;  // short tip, early termination of traversal
        nv = asg_arc_n(sg, v);
        av = asg_arc_a(sg, v);
        if ((v>>1)==(av[0].v>>1)) return 0;  // circle

        cnt++;
        // acc += (uint32_t)sg->arc[v].ul;
        // acc += (uint32_t) sg->idx[]

        for (i=0; i<nv; i++){
            if (!av[i].del) {
                v = av[0].v; 
                acc += (uint32_t) av[0].ul;
                break;
                }
        }
        if (hamt_asgarc_util_countPre(sg, v, 0, 0)!=1) return 0;  // backward bifur 
    }
    return 1;
}


int hamt_asgarc_util_get_the_one_target(asg_t *g, uint32_t v, uint32_t *w, int include_del_seq, int include_del_arc){
    // v as exactly one target (sucessor), given the `include_del_seq` and `include_del_arc` (NOT tested in this function)
    // get that vertex's index in av (NOT vertex ID)
    int nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    for (i=0; i<(int)nv; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        if (w){
            *w = av[i].v;
        }
        return i;
    }
    // fprintf(stderr, "[E::%s] unexpected.\n", __func__); exit(1);
    return -1;
}

int hamt_asgarc_util_walk_furthest_strict_single_path(asg_t *g, uint32_t v0, uint32_t *vn, int *nb_nodes, int limit){
    // FUNC
    //     walk the strict single path starting from v0
    //     store the target node ID in vn (if v0 has multiple targets, vn will just be v0)
    // RETURN
    //      1 if normal
    //      0 if reached limit
    //      -1 if the single path leads to cirling back to v0 (IN EITHER DIRECTION)
    uint32_t v=v0;
    asg_arc_t *av;
    *vn = v0;
    int cnt = 0;

    while (1){
        if (hamt_asgarc_util_countSuc(g, v, 0, 0)!=1){
            break;
        }
        hamt_asgarc_util_get_the_one_target(g, v, &v, 0, 0);
        if (hamt_asgarc_util_countPre(g, v, 0, 0)>0){
            break;
        }
        if ((v>>1)==(v0>>1)){
            return 0;
        }
        *vn = v;
        cnt++;
        if (cnt>=limit){
            break;
        }

    }
    if (nb_nodes){
        *nb_nodes = cnt;
    }
    return 1;    

}

int hamt_asgarc_util_countSinglePath(asg_t *g, uint32_t v0, int include_del_seq, int include_del_arc, 
                                     int max_arcs, int max_length, asg64_v *a, int *traversal_dist){
    // count length of single path; 
    // return 1 if encounters non-single edge within the given threshold
    // return 0 if longer than the given threshold
    // return -1 if loop
    // return -2 if meets end of the unitig
    //     for status code, use ha's *detect_single_path* functions.
    // (traverse in the direction of v. Feed with v^1 if looking for the opposite direction.)
    assert( !((max_arcs>0) && (max_length>0)) ); // not both
    assert( ((max_arcs>0) || (max_length>0)) ); // at least one
    int nb_fwd, nb_bwd;
    uint32_t v, w;
    int wid;
    asg_arc_t *av;
    v = v0;
    *traversal_dist = 0;  // init

    while (1){
        nb_fwd = hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc);
        if (nb_fwd==0) {return -2;}  // reaches the end of this tip
        if (nb_fwd>1) {break;}  // the end, or more than 1 targets
        wid = hamt_asgarc_util_get_the_one_target(g, v, &w, include_del_seq, include_del_arc);
        av = asg_arc_a(g, v);
        if (wid==-1){
            fprintf(stderr, "waht?\n");fflush(stderr);
            exit(1);
        }
        if ((w>>1) == (v0>>1)){
            return -1;  // loop
        }

        nb_bwd =  hamt_asgarc_util_countPre(g, w, include_del_seq, include_del_arc);
        assert(nb_bwd!=0);
        if (nb_bwd>1) {break;}  // w has more than 1 predecessor

        v = w;
        if (max_arcs>0) {
            *traversal_dist = *traversal_dist +1;
            if ((*traversal_dist)>max_arcs) {return 0;}
        }else {
            *traversal_dist += (int) ((uint32_t)av[wid].ul);
            if ((*traversal_dist)>max_length) {return 0;}
        }

        if (a){
            kv_push(uint64_t, *a, (uint64_t)w);
        }
    }
    return 1;  // encountered non-single edge

}

int hamt_asgarc_util_countNoneTipSuc(asg_t *g, uint32_t v){
    if (asg_arc_n(g, v)==0){
        return 0;
    }
    uint32_t w, nv;
    asg_arc_t *av;
    int cnt = 0;
    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        w = av[i].v;
        if (asg_arc_n(g, w)==0 && asg_arc_n(g, w^1)==1){
            continue;
        }
        cnt++;
    }
    return cnt;
}

int hamt_asgarc_util_isTip(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc){
    // FUNC
    //     check if v is a tip (or orphan), in either direction
    // RETURN
    //     0 if no
    //     1 if yes
    if (hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc)==0 || 
        hamt_asgarc_util_countPre(g, v, include_del_seq, include_del_arc)==0){
        return 1;
    }
    return 0;
}

int hamt_asgarc_util_checkSimpleBubble(asg_t *g, uint32_t v0, uint32_t *dest){
    // FUNC
    //     check if v0 is the start of a simple bubble (i.e. a 1-vertex bubble)
    // RETURN
    //     0 if no
    //     1 if yes
    int verbose = 0;
    if (hamt_asgarc_util_countSuc(g, v0, 0, 0)!=2){
        if (verbose){fprintf(stderr, "[checkBubble] v-suc\n");}
        return 0;
    }
    uint32_t w[2], nv, u, u1, u2;  // x is the next-next vertex (if exists)
    asg_arc_t *av = asg_arc_a(g, v0);
    int idx = 0;
    int idx_target;
    
    for (uint32_t i=0; i<asg_arc_n(g, v0); i++){
        if (av[i].del) {continue;}
        if (g->seq[av[i].v>>1].del) {continue;}
        w[idx] = av[i].v;
        idx++;
    }
    assert(idx==2);
    if (verbose){
        fprintf(stderr, "  (w1 utg%.6d w2 utg%.6d)\n", (int)(w[0]>>1)+1, (int)(w[1]>>1)+1);
    }

    // check if the next-next vertex is the end of a bubble
    if (hamt_asgarc_util_countSuc(g, w[0], 0, 0)!=1 || hamt_asgarc_util_countSuc(g, w[1], 0, 0)!=1){  // need to end in only one vertex
        if (verbose){fprintf(stderr, "[checkBubble] w-suc\n");}
        return 0;
    }
    if (hamt_asgarc_util_get_the_one_target(g, w[0], &u1, 0, 0)==-1){fprintf(stderr, "[%s] ??????", __func__);return 0;}
    if (hamt_asgarc_util_get_the_one_target(g, w[1], &u2, 0, 0)==-1){fprintf(stderr, "[%s] ??????", __func__);return 0;}
    if (u1!=u2){  // and it's the same vertex
        if (verbose){fprintf(stderr, "[checkBubble] u\n");}
        return 0;
    }
    if (hamt_asgarc_util_countPre(g, u1, 0, 0)!=2){  // and the "end of the simple bubble" must not have other incoming sources
        if (verbose){fprintf(stderr, "[checkBubble] u-pre (u is utg%.6d , pre: %d)\n", (int)(u1>>1)+1, hamt_asgarc_util_countPre(g, u1, 0, 0));}
        return 0;
    }
    // check backward links of the bubble edges
    if (hamt_asgarc_util_countPre(g, w[0], 0, 0)!=1 || hamt_asgarc_util_countPre(g, w[1], 0, 0)!=1){
        if (verbose){fprintf(stderr, "[checkBubble] w-suc\n");}
        return 0;
    }
    if (dest){
        *dest = u1;
    }
    return 1;
}

int hamt_asgarc_util_checkSimpleBubble_edge(asg_t *g, uint32_t v0){
    // check if v0 is an edge of a simple bubble (i.e. 1-vertex bubble)
    //     (doesn't allow tip, but could've; modify if needed)
    // RETURN
    //     0 if no
    //     1 if yes

    uint32_t s, e, e2, sib;  // start, end, sibling
    uint32_t nv;
    asg_arc_t *av;

    // check simple edge
    if (hamt_asgarc_util_countPre(g, v0, 0, 0)!=1 || hamt_asgarc_util_countSuc(g, v0, 0, 0)!=1){
        return 0;
    }
    // check both ends of the persumed bubble
    if (hamt_asgarc_util_get_the_one_target(g, v0, &e, 0, 0)<0){
        fprintf(stderr, "[ERROR::%s] e\n", __func__);fflush(stderr);
        exit(1);
    }
    if (hamt_asgarc_util_get_the_one_target(g, v0^1, &s, 0, 0)<0){
        fprintf(stderr, "[ERROR::%s] s\n", __func__);fflush(stderr);
        exit(1);
    }
    s = s^1;
    if (hamt_asgarc_util_countPre(g, e, 0, 0)!=2 || hamt_asgarc_util_countSuc(g, s, 0, 0)!=2){
        return 0;
    }
    // check the sibling
    nv = asg_arc_n(g, s);
    av = asg_arc_a(g, s);
    uint32_t i;
    for (i=0; i<nv; i++){
        if (av[i].v!=v0){
            sib = av[i].v;
            break;
        }
    }
    assert(i!=nv);
    if (hamt_asgarc_util_countPre(g, sib, 0, 0)!=1 || hamt_asgarc_util_countSuc(g, sib, 0, 0)!=1){
        return 0;
    }
    if (hamt_asgarc_util_get_the_one_target(g, sib, &e2, 0, 0)<0){
        fprintf(stderr, "[ERROR::%s] e2\n", __func__);fflush(stderr);
        exit(1);
    }
    if (e2==e){
        return 1;
    }
    return 0;
}

int hamt_asgarc_util_checkSimpleBubble_multiEddge(asg_t *sg, uint32_t v0, uint32_t *u0){
    // check if v0 is the start of a multi-edge simple bubble (including 2-edge simple bubble)
    // RETURN
    //      0 if no
    //      1 if yes (also right the end vertex to *u0)
    uint32_t w, u_tmp, u, nv, nw;
    int idx = 0;
    asg_arc_t *av, *aw;
    if (hamt_asgarc_util_countNoneTipSuc(sg, v0)<2){
        return 0;
    }
    nv = asg_arc_n(sg, v0);
    av = asg_arc_a(sg, v0);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (hamt_asgarc_util_isTip(sg, av[i].v, 0, 0)){continue;}
        w = av[i].v;
        if (hamt_asgarc_util_countPre(sg, w, 0, 0)!=1 || hamt_asgarc_util_countNoneTipSuc(sg, w)!=1){
            return 0;
        }
        // check each edge's target
        aw = asg_arc_a(sg, w);
        for (uint32_t i2=0; i2<asg_arc_n(sg, w); i2++){  // find the non-tip target
            if (aw[i2].del){continue;}
            if (hamt_asgarc_util_isTip(sg, aw[i2].v, 0, 0)){continue;}
            u_tmp = aw[i2].v;
            if (idx==0){
                u = u_tmp;
            }else{
                if (u_tmp!=u){
                    return 0;
                }
            }
            idx++;
            break;
        }
    }
    if (idx!=hamt_asgarc_util_countNoneTipSuc(sg, u^1)){
        return 0;
    }
    *u0 = u;
    return 1;
}

int hamt_asgarc_util_checkSimpleInvertBubble(asg_t *sg, uint32_t v0, uint32_t *w0, uint32_t *u0){
    /*
        v0->w1
        v0->u1
        u1->....
        w1->u1^1  
        we want to drop v0->w1 and w1->u1^1
    */
   // only treat the most simple case, any branching and it'll be ignored
   // RETURN
   //     1 if v0 is the startpoint of such structure
   //     0 otherwise
    uint32_t uw[2], x;
    int nuw[4];
    uint32_t nv, idx=0;
    asg_arc_t *av;
    if (hamt_asgarc_util_countSuc(sg, v0, 0, 0)!=2){
        return 0;
    }
    nv = asg_arc_n(sg, v0);
    av = asg_arc_a(sg, v0);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){
            continue;
        }
        uw[idx] = av[i].v;
        idx++;
    }
    assert(idx==2);
    
    // check topo
    nuw[0] = hamt_asgarc_util_countSuc(sg, uw[0], 0, 0);
    nuw[1] = hamt_asgarc_util_countSuc(sg, uw[1], 0, 0);
    nuw[2] = hamt_asgarc_util_countPre(sg, uw[0], 0, 0);
    nuw[3] = hamt_asgarc_util_countPre(sg, uw[1], 0, 0);
    if (nuw[0]==0 || nuw[1]==0){  // no target
        return 0;
    }
    if (nuw[0]>1 && nuw[1]>1){  // only one of them can have more than 1 target
        return 0;
    }
    if (nuw[0]==1 && nuw[1]==1){  // no inversion possible
        return 0;
    }
    if (nuw[2]>1 || nuw[3]>1){  // branching in backward direction
        return 0;
    }
    // check topo: the inversion
    if (nuw[0]==1){
        if (hamt_asgarc_util_get_the_one_target(sg, uw[0], &x, 0, 0)==-1){
            fprintf(stderr, "[W::%s] failed to get target when supposed to.\n", __func__);
            return 0;
        }
        if (x!=(uw[1]^1)){
            return 0;
        }else{
            *w0 = uw[0];
            *u0 = uw[1];
            return 1;
        }
    }else{
        assert(nuw[1]==1);
        if (hamt_asgarc_util_get_the_one_target(sg, uw[1], &x, 0, 0)==-1){
            fprintf(stderr, "[W::%s] failed to get target when supposed to.\n", __func__);
            return 0;
        }
        if (x!=(uw[0]^1)){
            return 0;
        }else{
            *w0 = uw[1];
            *u0 = uw[0];
            return 1;
        }
    }

}


int hamt_asgarc_util_checkDontCovCut_givenArc(asg_t *g, asg_arc_t *arc){
    // NOTE
    //       intented to be used by SCCcovCut
    // FUNC
    //    check if the arc is in a (slightly relaxed) simple bubble, or end of a tip
    //      (allow more than 2 bubble edges, but no other more complex structures)
    // RETURN
    //    0 if don't cut
    //    1 if yes consider cutting
    uint32_t v = arc->ul>>32, w = arc->v, u, ww;
    uint32_t nv = hamt_asgarc_util_countSuc(g, v, 0, 0);
    uint32_t nw = hamt_asgarc_util_countSuc(g, w, 0, 0);
    uint32_t nw_ = hamt_asgarc_util_countPre(g, w, 0, 0);
    asg_arc_t *av, *aw_;
    if (nv>1 && nw>1){  // doesn't look like a bubble
        return 1;
    }
    if (nw==0 || hamt_asgarc_util_countPre(g, v, 0, 0)==0){  // is a tip end, dont cut
        return 0;
    }
    // might be a bubble
    if (nv>1){  // v be the start
        assert(nw==1);
        if (hamt_asgarc_util_countPre(g, w, 0, 0)>1){
            return 1;
        }
        u = asg_arc_a(g, w)[0].v;
        av = asg_arc_a(g, v);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            ww = av[i].v;
            if (hamt_asgarc_util_countSuc(g, ww, 0, 0)>1 || hamt_asgarc_util_countPre(g, ww, 0, 0)>1){  // branching in either direction inside the persumed bubble
                return 1;
            }
            uint32_t tmpv;
            if (hamt_asgarc_util_get_the_one_target(g, ww, &tmpv, 0, 0)>0){
                if (tmpv!=u){  // persumed bubble edge goes to another place
                    return 1;
                }
            }
        }
        return 0;
    } else {  // w be the end
        assert(nv==1);
        if (hamt_asgarc_util_countPre(g, v, 0, 0)>1){
            return 1;
        }
        u = asg_arc_a(g, v^1)[0].v;
        aw_ = asg_arc_a(g, w^1);
        for (uint32_t i=0; i<nw_; i++){
            ww = aw_[i].v;
            if (hamt_asgarc_util_countSuc(g, ww, 0, 0)>1 || hamt_asgarc_util_countPre(g, ww, 0, 0)>1){
                return 1;
            }
            uint32_t tmpv;
            if (hamt_asgarc_util_get_the_one_target(g, ww, &tmpv, 0, 0)>0){
                if (tmpv!=u){  // persumed bubble edge goes to another place
                    return 1;
                }
            }
        }
        return 0;
    }

}

int hamt_asgarc_util_countSinglePath_allowSimpleBubble(asg_t *g, uint32_t v0, int max_arcs){
    // REQUIRE
    //     graph is clean of loops
    // FUNC
    //     count length of single path, allowing 1-vertex bubbles
    // RETURN
    //     return 1 if terminated because of branching before max_arcs
    //     return 0 if exceeded max_arcs
    //     return -1 if terminated because of end of subgraph
    assert(max_arcs>0);
    uint32_t v = v0, w, w1, w2;
    int nb_fwd, nb_bwd;  // suf/pre counts
    int dist = 0;
    while (1){
        nb_fwd = hamt_asgarc_util_countSuc(g, v, 0, 0);
        if (nb_fwd==0) {  // reached the end of this single path
            return -1;
        } else if (nb_fwd==1){
            hamt_asgarc_util_get_the_one_target(g, v, &w, 0, 0);
            nb_bwd = hamt_asgarc_util_countPre(g, w, 0, 0);
            if (nb_bwd!=1){
                return 1;  // the next target has other incoming sources, terminate
            }else{  // step
                v = w;
                dist++;
            }
        } else if (nb_fwd==2){
            int ret = hamt_asgarc_util_checkSimpleBubble(g, v, &w);
            if (ret){
                dist = dist+2; // because bubble
                v = w;
            }else{
                return 1; // not a simple bubble, for whatever reason
            }
        }else{
            return 1;
        }
        if (dist>max_arcs){
            return 0;
        }
    }
    fprintf(stderr, "%s: unexpected\n", __func__); fflush(stderr);
    exit(1);
}

int hamt_ugasgarc_util_BFS_mark(ma_ug_t *ug, uint32_t v0, uint8_t *color, queue32_t *q,
                                int threshold_nread, int threshold_bp){
    // helper func, do a BFS within the given range(s) and mark touched vertices in `color`
    // (if both thresholds are given, require exceed them both to terminate search)
    // (set threshold to -1 to disable)
    // the caller is responsible for init color and q
    // RETURN
    //     1 if terminated by emptying the queue (early termination)
    //     0 if terminated by reaching the limit
    asg_t *auxsg = ug->g;
    uint32_t v, nv, w;
    asg_arc_t *av;
    int cnt_nread=0, cnt_bp=0;
    int tmp;
    queue32_enqueue(q, v0);
    color[v0] = 1;
    int is_early_termination = 1;
    while(queue32_dequeue(q, &v)){
        // fprintf(stderr, "dequeue, %d, queue size is %d, start %d, end %d\n", (int)(v>>1)+1, q->n, (int)q->i_head, (int)q->i_tail);
        if (color[v]==2){
            continue;
        }
        nv = asg_arc_n(auxsg, v);
        av = asg_arc_a(auxsg, v);
        tmp = 0;  // used to get the length of the longest child utg 
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            w = av[i].v;
            if (color[w]==0){
                queue32_enqueue(q, w);
                // fprintf(stderr, "  push, %d\n", (int)(w>>1)+1);
                color[w] = 1;
                tmp = tmp<ug->u.a[w>>1].len? ug->u.a[w>>1].len : tmp;
            }
            // this isn't the length of the "straight" path of course, 
            // but will serve better if there's a hairball
            cnt_nread++;
            cnt_bp += tmp;  
        }
        color[v] = 2;

        if (threshold_bp<0 && threshold_nread<0){
            continue;
        }
        if (threshold_bp>0 && threshold_nread<0 &&cnt_bp>threshold_bp){
            is_early_termination = 0;
            break;
        }
        if (threshold_nread>0 && threshold_bp<0 && cnt_bp>threshold_nread){
            is_early_termination = 0;
            break;
        }
        if (threshold_nread>0 && threshold_bp>0 && cnt_nread>threshold_nread && cnt_bp>threshold_bp){
            is_early_termination = 0;
            break;
        }
    }
    return is_early_termination;
}


int hamt_asgarc_util_BFS_markSubgraph(asg_t *sg, int *subg_labels, uint8_t *color, queue32_t *q){
    // helper func, do a BFS within the given range(s) and mark touched vertices in `color`
    // (if both thresholds are given, require exceed them both to terminate search)
    // (set threshold to -1 to disable)
    // the caller is responsible for init color and q
    int verbose = 0;

    uint32_t v0, v, nv, w;
    asg_arc_t *av;
    int cnt_nread, cnt_bp;
    int tmp;
    int label = 0;

    for (v0=0; v0<sg->n_seq*2; v0++){
        if (color[v0]){
            continue;
        }
        tmp = 0;
        queue32_enqueue(q, v0);
        color[v0] = 1;
        tmp++;
        subg_labels[v0] = label;
        while(queue32_dequeue(q, &v)){
            if (color[v]==2){
                continue;
            }
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            for (uint32_t i=0; i<nv; i++){
                // if (av[i].del){continue;}
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(q, w);
                    color[w] = 1;
                    tmp++;
                    subg_labels[w] = label;
                }
            }
            color[v] = 2;
        }
        label++;
        if (verbose){
            fprintf(stderr, "[debug::%s] subgraph #%d, %d reads.\n", __func__, label, tmp);
        }
    }
    return label;
}

int hamt_asgutil_calc_singlepath_cov(asg_t *g, uint32_t v){
    // FUNC
    //     given v, find the kmer-based coverage of the longest single path IN THE REVERSE DIRECTION 
    //     (because it's intend to be sused with hamt_asgutil_detect_dangling_circle)

    uint32_t w = v^1;
    int wid;
    float cov=0;
    int n=0;
    while (1){
        // log the current vertex
        n++;
        cov+=R_INF.median[w>>1];

        // step to the next one
        if (hamt_asgarc_util_countSuc(g, w, 0, 0)!=1){
            break;
        }
        wid = hamt_asgarc_util_get_the_one_target(g, w, &w, 0, 0);
        // w = asg_arc_a(g, w)[wid].v;
        if (hamt_asgarc_util_countPre(g, w, 0, 0)!=1){
            break;
        }
    }
    cov = cov/n;
    return (int) (cov+0.499);
}

int hamt_asgutil_is_tigSafeToErode(asg_t *g, uint32_t v){
    // REQUIRE
    //     v is a tip, in the direction of: has no predecessor, has target(s)
    // FUNC
    //     we don't want to discard small unitigs all together when eroding the tips,
    //     so his function checks whether a given vertex is in a relative small single path (allow bubbles) 
    // NOTE
    //    doesn't protect branched small subgraphs - they'll only be left with the "backbone" (which is randomly decided)
    // RETURN
    //    1 if yes
    //    0 if no
    int has_branching;
    has_branching = hamt_asgarc_util_countSinglePath_allowSimpleBubble(g, v, 10);  // the number of vertices, not bp
    if (has_branching==1){
        return 1;  // this tip can be erode
    }else{
        return 0;  // don't touch it
    }

}


// int hamt_asgarc_both_in_coherent_circle(asg_t *sg, uint32_t v, uint32_t w){
//     // NOTE
//     //     not a good way to do this...
//     //     this is a hot patch thing for SCC-coverage-cut since SCC labeling has been weirdly not right
//     //     if SCC improves, maybe ditch this func
//     // FUNC
//     //     v->w forms an arc. This function tests if there's a coherent path
//     //     (i.e. no incoherent edge-vertex-edge spot in terms of orientation)
//     //     starts from w and goes to v, i.e. v and w is in a circle
//     // TODO
//     //     if keeping this func: alloc buffers in the caller.

//     // simplest cases
//     uint32_t nv;
//     asg_arc_t *av;
//     nv = asg_arc_n(sg, w);
//     av = asg_arc_a(sg, w);
//     for (int i=0; i<nv; i++){
//         if (av[i].v==v){
//             return 1;
//         }
//     }

//     // traversal buffer
//     int8_t *orientation = (int8_t*)malloc(sg->n_seq*2*1);
//     memset(orientation, -1, sg->n_seq*2);
//     uint8_t *color = (uint8_t*)calloc(sg->n_seq*2, 1);
//     stacku32_t stack;
//     stacku32_init(&stack);

//     // 


//     free(orientation);
//     free(color);
//     stacku32_destroy(&stack);
    
// }

///////////////////////////////////////////////////////////////////////////////////////
//                  higher-level routines (assembly graph)                           //
///////////////////////////////////////////////////////////////////////////////////////


asg_t *hamt_asggraph_util_gen_transpose(asg_t *g0){
    // note that not everything is copied; 
    // this function is only intended to be used by DFS/SCC
    asg_t *g = asg_init();
    
    g->m_arc = g0->m_arc;
    g->n_arc = g0->n_arc;
    g->is_srt = 0;
    g->m_seq = g0->m_seq;
    g->n_seq = g0->n_seq;
    g->is_symm = 1;

    g->arc = (asg_arc_t*)malloc(sizeof(asg_arc_t) * g0->n_arc);
    asg_arc_t *a, *a0;
    for (uint32_t i=0; i<g->n_arc; i++){
        a = &g->arc[i];
        a0 = &g0->arc[i];
        a->ul = ((uint64_t)a0->v)<<32 | (uint64_t)((uint32_t)a0->ul);
        a->v = a0->ul>>32;
    }
    
    asg_arc_sort(g);
    g->is_srt = 1;
    asg_arc_index(g);

    return g;

}


void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length){
    // rationale: tips could hinder simple bubble detection, 
    //            and hamt doesn't need to be too careful with both of them (i think) (Sep 17)  
    double startTime = Get_T();
    uint32_t n_vtx = g->n_seq*2, v;
    uint32_t n_del;
    int pre, suc, l, ret;

    asg64_v a = {0,0,0};
    
    for (int i_round=0; i_round<3; i_round++)
    {
        n_del = 0;

        // drop simple circles
        n_del += asg_arc_del_simple_circle_untig(NULL, NULL, g, 100, 0);

        // pop bubbles
        // (comment in ha was "remove isoloated single read")
        n_del += asg_arc_del_single_node_directly(g, asm_opt.max_short_tip, sources);

        // tip
        n_del += asg_cut_tip(g, asm_opt.max_short_tip);

        // check
        if (n_del==0)
            {break;}
        else{
            asg_cleanup(g);
            if (VERBOSE>=1){
                fprintf(stderr, "[M::%s] removed %d locations\n", __func__, n_del);
                fprintf(stderr, "[M::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
            }
        }
        
        n_del = 0;
        n_del = hamt_sg_pop_simpleInvertBubble(g);
        if (VERBOSE>=1){
            fprintf(stderr, "[M::%s] removed %d invert bubble on string graph\n", __func__, n_del);
            fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
        }
    }
    free(a.a);

}






///////////////////////////////////////////////////////////////////////////////////////
//                  higher-level routines (unitig graph)                             //
///////////////////////////////////////////////////////////////////////////////////////

// TODO proper bi-directional SCC

uint32_t* hamt_ugarc_util_v2vu_nodirection(ma_ug_t *ug, asg_t *sg){
    uint32_t *v2vu = (uint32_t*)calloc(sg->n_seq, sizeof(uint32_t));  assert(v2vu);  // if you segfault here: did you forgot to >>1
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        v2vu[ug->u.a[vu].start>>1] = vu;
        v2vu[ug->u.a[vu].end>>1] = vu;
    }
    return v2vu;
}

uint32_t* hamt_ugarc_util_v2vu_yesdirection(ma_ug_t *ug, asg_t *sg){
    uint32_t *v2vu = (uint32_t*)calloc(sg->n_seq*2, sizeof(uint32_t));  assert(v2vu);
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        v2vu[ug->u.a[vu].start] = vu<<1 | 0;
        v2vu[ug->u.a[vu].end] = vu<<1 | 1;
    }
    return v2vu;
}

int hamt_asgarc_detect_circleDFS(asg_t *sg, uint32_t v0, uint32_t w0, int allow_v0){ // vu0 has direction; v2vu has directions
    // DFS on sg (could be ug's auxsg), test if vu is the start of a circle
    // (note "start of a circle" means v0 is the "stem"; v0 shall not be part of the circle)
    // return 1 if yes, 0 otherwise

    uint32_t v=w0, nv, w, wu;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(sg->n_seq*2, 1);
    stacku32_t stack;
    stacku32_init(&stack);
    stacku32_push(&stack, v);
    color[v] = 1;

    int ret = 0;

    // fprintf(stderr, ">checking utg%.6d\n", (int)(vu>>1)+1);
    while (1){  // DFS main loop
        if (stacku32_pop(&stack, &v)){
            // fprintf(stderr, "    pop utg%.6d (dir %d)\n", (int)(vu>>1)+1, (int)(vu&1));
            if (color[v]==2){
                continue;
            }
            // fprintf(stderr, "        handle is %.*s\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1));
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            if (hamt_asgarc_util_countSuc(sg, v, 0, 0)==0){
                color[v] = 2;
                continue;
            }
            int is_exhausted = 1;
            for (int i=0; i<(int)nv; i++){
                if (av[i].del){continue;}
                w = av[i].v;

                if (w==w0){  // found the circle
                    ret = 1;  // don't exit, we need to empty the stack to make sure the handle tig isn't in
                }
                if (w==v0){  // the circle shall not include the handle tig
                    if (!allow_v0){
                        // fprintf(stderr, "HANDLE IN CIRCLE.\n");
                        free(color);
                        stacku32_destroy(&stack);
                        return 0;
                    }else{
                        ;
                    }
                }

                if (color[w]==0){
                    if (is_exhausted){
                        is_exhausted = 0;
                        stacku32_push(&stack, v); // put the current node back
                        // fprintf(stderr, "    put back utg%.6d\n", (int)(vu>>1)+1);
                    }
                    color[w] = 1;
                    stacku32_push(&stack, w);
                    // fprintf(stderr, "    pushed utg%.6d (dir %d)\n", (int)(wu>>1)+1, (int)(wu&1));
                    // fprintf(stderr, "        handle is %.*s\n", (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1));
                }
            }
            if (is_exhausted){
                color[v] = 2;
            }
        }else{
            break;
        }
    }
    free(color);
    stacku32_destroy(&stack);
    return ret;
}

int hamt_asgarc_detect_circleDFS2(asg_t *sg, uint32_t v0, uint32_t w0, uint32_t *v_exclude){
    // FUNC
    //     check if v0 and w0 are in the same circle
    //     if v_exclude is not NULL, it points to ONE uint32_t, and it shall not appear during the traversal
    // RETURN
    //     1 if yes, 0 otherwise
    //     -1 if there might be the inversion case
    // NOTE
    //     might need 2 DFS.
    stacku32_t stack;
    stacku32_init(&stack);
    uint8_t *color = (uint8_t*)calloc(sg->n_seq*2, 1);
    int has_found_w0 = 0, has_found_v0 = 0, is_weird = 0, is_exhausted;
    int is_initial_pop = 1;
    uint32_t v, w, nv;
    asg_arc_t *av;

    // 1st direction
    stacku32_push(&stack, v0);
    color[v0] = 1;
    while (1){
        if (stacku32_pop(&stack, &v)){
            if (color[v]==2){continue;}
            if (!is_initial_pop){
                if (v==w0){ has_found_w0 = 1; }
                if (v==v0){ has_found_v0 = 1;}
                if (v==(v0^1) || v==(w0^1)){ is_weird++; }
            }else{
                is_initial_pop = 0;
            }
            if (v_exclude && (((*v_exclude)>>1)==(v>>1)) ){
                stacku32_destroy(&stack);
                free(color);
                return 0;
            }
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            is_exhausted = 1;
            for (uint32_t i=0; i<nv; i++){
                if (color[av[i].v]==0){
                    is_exhausted = 0;
                    break;
                }
            }
            if (!is_exhausted){
                stacku32_push(&stack, v);
                for (uint32_t i=0; i<nv; i++){
                    w = av[i].v;
                    if (color[w]==0){
                        stacku32_push(&stack, w);
                        color[w] = 1;
                    }
                }
            }else{
                color[v] = 2;
            }
        }else{
            break;
        }

    }
    if (!has_found_w0 || !has_found_v0){
        stacku32_destroy(&stack);
        free(color);
        return 0;
    }
    if (is_weird>=2){
        stacku32_destroy(&stack);
        free(color);
        return -1;
    }

    // reset
    stacku32_reset(&stack);
    memset(color, 0, sg->n_seq*2*1);
    stacku32_push(&stack, v0^1);
    color[v0^1] = 1;
    is_initial_pop = 1;
    has_found_v0 = 0;
    has_found_w0 = 0;

    // 2nd direction
    while (1){
        if (stacku32_pop(&stack, &v)){
            if (color[v]==2){continue;}
            if (!is_initial_pop){
                if (v==(w0^1)){ has_found_w0 = 1; }
                if (v==(v0^1)){ has_found_v0 = 1;}
            }else{
                is_initial_pop = 0;
            }
            if (v_exclude && (((*v_exclude)>>1)==(v>>1)) ){
                stacku32_destroy(&stack);
                free(color);
                return 0;
            }
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            is_exhausted = 1;
            for (uint32_t i=0; i<nv; i++){
                if (color[av[i].v]==0){
                    is_exhausted = 0;
                    break;
                }
            }
            if (!is_exhausted){
                stacku32_push(&stack, v);
                for (uint32_t i=0; i<nv; i++){
                    w = av[i].v;
                    if (color[w]==0){
                        stacku32_push(&stack, w);
                        color[w] = 1;
                    }
                }
            }else{
                color[v] = 2;
            }
        }else{
            break;
        }
    }
    stacku32_destroy(&stack);
    free(color);
    if (has_found_w0 && has_found_v0){
        return 1;
    }else{
        return 0;
    }    
}

void hamt_asgarc_ugCovCutDFSCircle(asg_t *sg, ma_ug_t *ug)
{
    int verbose = 1;  // debug
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t v, w, nv;
    // uint32_t cov1, cov2, covtmp, cov3, covmin,;
    uint32_t cov[3];
    uint32_t covtmp;
    asg_arc_t *av;
    // uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));  // for utg coverage

    int nb_cut =0;

    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del){
            continue;
        }
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (hamt_asgarc_util_countPre(auxsg, w, 0, 0)<2){
            continue;
        }

        cov[2] = 0;
        av = asg_arc_a(auxsg, w^1);
        for (uint32_t i=0; i<asg_arc_n(auxsg, w^1); i++){  // collect the largest coverage of w's predecessor (except v which will be accounted by cov[0])
            if ((av[i].v>>1)==(v>>1)){  // ignore v
                continue;
            }
            if (av[i].del){continue;}
            // covtmp = get_ug_coverage(&ug->u.a[av[i].v>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
            covtmp = ug->utg_coverage[av[i].v>>1];
            if (covtmp>cov[2]){
                cov[2] = covtmp;
            }
        }
        // cov[0] = get_ug_coverage(&ug->u.a[v>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
        cov[0] = ug->utg_coverage[v>>1];
        // cov[1] = get_ug_coverage(&ug->u.a[w>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
        cov[1] = ug->utg_coverage[w>>1];
        radix_sort_ovhamt32(cov, cov+3);

        // check coverage
        if (!hamt_check_covDiff(cov[0], cov[2])){
            fprintf(stderr, "[debug::%s] failed covDiff utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
            continue;
        }

        // detect circles and cut
        if (hamt_asgarc_detect_circleDFS(auxsg, v, w, 0)){  // note: v and w is ug with direction
            hamt_ug_arc_del(sg, ug, v, w, 1);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
            }
        }

    }
    
    if (nb_cut){
        asg_cleanup(sg); asg_symm(sg);
        asg_cleanup(auxsg); asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d links.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}


void hamt_asgarc_ugCovCutDFSCircle_aggressive(asg_t *sg, ma_ug_t *ug){
    // FUNC
    //      when performing DFS-based circle detection,
    //      dont require the handle to be not included in the circle
    //      criteria
    //         1) there exists a circling path, without touching the handle
    //         2) the circling path involving the handle shall be much longer than the one in 1)
    //         3) more relaxed pre/suc count requirements and coverage diff requirements
    // NOTE
    //     experimental, only use it after all more carefule steps have been done
    int verbose = 1;  // debug
    double startTime = Get_T();
    int nb_cut = 0;

    asg_t *auxsg = ug->g;
    uint32_t v, w, nv, nw;
    asg_arc_t *av, *aw;
    uint32_t covs[20], cov_min, cov_max, cov_v, cov_tmp;  // TODO: tho it's probably safe enough

    for (v=0; v<auxsg->n_seq*2; v++){
        // skip tips and simple bubble starts
        if (hamt_asgarc_util_countNoneTipSuc(auxsg, v)==0){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleBubble(auxsg, v, NULL)){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleBubble_edge(auxsg, v)){
            continue;
        }

        av = asg_arc_a(auxsg, v);
        nv = asg_arc_n(auxsg, v);
        for (int iw=0; iw<nv; iw++){
            if (av[iw].del){continue;}
            w = av[iw].v;
            if (hamt_asgarc_util_isTip(auxsg, w, 0, 0)){
                continue;
            }
            if (hamt_asgarc_util_countPre(auxsg, w, 0, 0)<2){
                continue;
            }
            // check coverage diff
            nw = asg_arc_n(auxsg, w^1);
            aw = asg_arc_a(auxsg, w^1);
            int idx = 0;
            for (int i=0; i<nw; i++){
                if (aw[i].del){continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (hamt_asgarc_util_isTip(auxsg, aw[i].v, 0, 0)){continue;}  // don't use tip's coverage
                covs[idx] = (uint32_t)ug->utg_coverage[aw[i].v>>1];
                idx++;
            }
            covs[idx++] = (uint32_t)ug->utg_coverage[w>>1];
            radix_sort_ovhamt32(covs, covs+idx);
            cov_min = covs[0];
            cov_max = covs[idx-1];
            cov_v = ug->utg_coverage[v>>1];
            int yes_cut = 0;
            if (cov_v<cov_min){
                if ((cov_max-cov_v)>10){yes_cut = 1;}
            }else if (cov_v<cov_max){
                if ((cov_max-cov_v)>10 || (cov_v-cov_min)>10){yes_cut = 1;}
            }else{
                if ((cov_v-cov_min)>10){yes_cut = 1;}
            }
            if (!yes_cut){
                if (verbose){
                    fprintf(stderr, "[debug::%s] failed cov diff, utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
                continue;
            }

            yes_cut = 0;
            // passed coverage diff check, detect circling paths
            // special case: 1-vertex circle
            for (int i=0; i<nw; i++){
                if (aw[i].del){continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (aw[i].v==(w^1)){
                    yes_cut = 1;
                    break;
                }
            }
            if (yes_cut){
                hamt_ug_arc_del(sg, ug, v, w, 1);
                nb_cut++;
                if (verbose){
                    fprintf(stderr, "[debug::%s] special cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
                continue;
            }
            // regular case
            for (int i=0; i<nw; i++){
                if (aw[i].del){continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (hamt_asgarc_detect_circleDFS2(auxsg, w^1, aw[i].v, &v)>0){
                    yes_cut = 1;
                    break;
                }
            }
            if (yes_cut){
                hamt_ug_arc_del(sg, ug, v, w, 1);
                nb_cut++;
                if (verbose){
                    fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
                continue;
            }else{
                if (verbose){
                    fprintf(stderr, "[debug::%s] DIDNOT cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
            }
        }
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    }

}

int hamt_asgarc_contain_only_simple_circle(asg_t *sg, asg_t *auxsg, int *SCC_marking, int SCC_label){
    // NOTE
    //     intented to by used by SCC-coverage-cut
    queue32_t q;
    queue32_init(&q);
    uint32_t v, w, w_, nv, v0;
    int seeded = 0, passed = 0;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    int suc_cnt, suc_cnt_;
    int terminate = 0;

    for (v=0; v<auxsg->n_seq*2; v++){
        if (SCC_marking[v]==SCC_label){
            if (hamt_asgarc_util_countSuc(auxsg, v, 0, 0)>0){
                seeded = 1;
                v0 = v;
                break;
            }
        }
    }
    if (!seeded){
        free(color);
        queue32_destroy(&q);
        return 0;
    }
    queue32_enqueue(&q, v);
    color[v] = 1;
    while (queue32_dequeue(&q, &v)){
        if (v==v0){
            passed = 1;  // dont terminate yet, we need to empty the whole queue
        }
        if (color[v]==2){
            continue;
        }
        nv = asg_arc_n(auxsg, v);
        av = asg_arc_a(auxsg, v);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            w = av[i].v;
            w_ = w^1;
            if (color[w]==0){
                // exclude cases with weird edges
                if (color[w^1]!=0){
                    terminate = 1;
                    break;
                }
                // check topo
                suc_cnt = hamt_asgarc_util_countNoneTipSuc(auxsg, w);
                suc_cnt_ = hamt_asgarc_util_countNoneTipSuc(auxsg, w_);
                // (check forward)
                if (suc_cnt==0){  // no none-tip successor
                    continue;
                }
                if (suc_cnt>2){  // not a simple circle
                    terminate = 1;
                    break;
                }
                if (suc_cnt==2){  // check if it's the start of a simple bubble
                    if (!hamt_asgarc_util_checkSimpleBubble(auxsg, w, NULL)){
                        terminate = 1;
                        break;
                    }
                }
                // check backward
                if (suc_cnt_>2){terminate = 1; break;}  // not a simple circle
                if (suc_cnt_==2){
                    if (!hamt_asgarc_util_checkSimpleBubble(auxsg, w^1, NULL)){
                        terminate = 1;
                        break;
                    }
                }

                // push
                color[w] = 1;
                queue32_enqueue(&q, w);
            }
        }
        color[v] = 2;

        if (terminate){
            queue32_destroy(&q);
            free(color);
            return 0;
        }
    }

    queue32_destroy(&q);
    free(color);
    return 1;
}


int hamt_asgarc_remove_inter_SCC(asg_t *sg, asg_t *auxsg, ma_ug_t *ug, int *SCC_marking, int SCC_label){
    // NOTE
    //     dont cut the tips
    uint32_t v, w, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (SCC_marking[v]!=SCC_label && SCC_marking[w]!=SCC_label){  // irrelevant
            continue;
        }
        if (SCC_marking[v]==SCC_label && SCC_marking[w]==SCC_label){  // inner arc
            continue;
        }
        // check tips
        if (SCC_marking[v]==SCC_label){
            if (hamt_asgarc_util_countSuc(auxsg, w, 0, 0)==0){
                continue;
            }
        }else{
            if (hamt_asgarc_util_countPre(auxsg, v, 0, 0)==0){
                continue;
            }
        }
        // cut on the string graph
        hamt_ug_arc_del(sg, ug, v, w, 1);
        // if ((v&1)==0){
        //     v = ug->u.a[v>>1].end^1;
        // }else{
        //     v = ug->u.a[v>>1].start^1;
        // }
        // if ((w&1)==0){
        //     w = ug->u.a[w>>1].start;
        // }else{
        //     w = ug->u.a[w>>1].end;
        // }
        // asg_arc_del(sg, v, w, 1);
        // asg_arc_del(sg, w^1, v^1, 1);
        nb_cut++;
    }
    return nb_cut;
}

int hamt_asgarc_rescueArcWithinSubgraph(asg_t *sg, ma_ug_t *ug, int *subg_labels, int label){
    // (experimental) recover previous deleted arcs within the subgraph 
    uint32_t v, w, vu, wu;
    int cnt = 0;
    if (!ug){
        for (uint32_t i=0; i<sg->n_arc; i++){
            if (!sg->arc[i].del){
                continue;
            }
            v = sg->arc[i].ul>>32;
            w = sg->arc[i].v;
            if (subg_labels[v]!=label || subg_labels[w]!=label){
                continue;
            }
            // if (hamt_asgarc_util_countSuc(sg, v, 0, 0)==0 && hamt_asgarc_util_countPre(sg, w, 0, 0)==0)
            {
                asg_arc_del(sg, v, w, 0);
                asg_arc_del(sg, w^1, v^1, 0);
                cnt++;
            }
        }
    }else{  // has unitig graph, rescue at ug level (modify both sg and auxsg)
        asg_t *auxsg = ug->g;
        for (uint32_t i=0; i<auxsg->n_arc; i++){
            if (!auxsg->arc[i].del){
                continue;
            }
            vu = auxsg->arc[i].ul>>32;
            wu = auxsg->arc[i].v;
            if (subg_labels[vu]!=label || subg_labels[wu]!=label){
                continue;
            }

            {
                hamt_ug_arc_del(sg, ug, vu, wu, 0);
                cnt++;
            }
        }

    }
    return cnt;
}

int hamt_asgarc_rescueArcWithinCurrentSubgraph(asg_t *sg, ma_ug_t *ug, uint32_t v0){
    // dont need to label all subgraphs, just search based on the given seed
    int verbose = 1;
    asg_t *auxsg = ug->g;
    int cnt = 0;

    // get the current subgraph
    uint32_t v, nv, w;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    queue32_t q;
    queue32_init(&q);
    uint32_t v02[2] = {v0, v0^1};
    for (int i=0; i<2; i++){
        v = v02[i];
        queue32_enqueue(&q, v);
        color[v] = 1;
        while(queue32_dequeue(&q, &v)){
            if (color[v]==2){
                continue;
            }
            nv = asg_arc_n(auxsg, v);
            av = asg_arc_a(auxsg, v);
            for (uint32_t i=0; i<nv; i++){
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(&q, w);
                    color[w] = 1;
                }
            }
            color[v] = 2;
        }
    }

    // rescue
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (!auxsg->arc[i].del){continue;}
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (!color[v] || !color[w]){continue;}
        hamt_ug_arc_del(sg, ug, v, w, 0);  // set arc to not deleted
        if (verbose){
            fprintf(stderr, "[debug::%s] rescued the arc between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
        }
        cnt++;
    }

    // clean up
    free(color);
    queue32_destroy(&q);

    return cnt;  // count of rescue (by spot)
}

/////////////////// bidirectional graph SCC (broken, dont use until proven otherwise)/////////////////////
/*
note on the SCC thing: (Oct 25 2020)
    Kazutoshi Ando et al. (1996) proved that bidirectional graph can be uniquely decomposed into SCCs (defined 
    similarly to the SCC in ordinarly directed graphs aka informally: there exist paths P1!=0, P2!=0 and P1+P2=0 
    that describes u->v and v->u, then u and v are strongly connected to each other, and this is an equivalence relation on the graph) 
    and this could be done in linear time by doing Tarjan's SCC algo on a (ordinary) directed graph that is derived 
    from the bidirectional graph. In hifiasm/miniasm's implementation, the derived auxiliary graph is asg_t *ug->g.

    However, under this definition, some intuitively "isolated" subgraphs would fall into the same SCC if there's 
    inconsistent nodes / inversions. For hamt heuristics, turns out this was not exactly the approach that was called for, and 
    bridge detection + DFS (although might be relatively expensive) might serve better.

    Code blocks have been removed.
*/
#if 0
// was the main handle
int *hamt_asg_SCC(asg_t *sg){  
    uint64_t *ts1 = (uint64_t*)calloc(sg->n_seq*2, sizeof(uint64_t)); assert(ts1);
    int *SCC_labels = (int*)calloc(sg->n_seq*2, sizeof(int));  assert(SCC_labels);

    hamt_asg_DFS(sg, ts1, NULL, 0, NULL);
    radix_sort_ovhamt64(ts1, ts1+sg->n_seq*2);
    hamt_asg_DFS(sg, NULL, ts1, 1, SCC_labels);

    free(ts1);
    return SCC_labels;
}
#endif
/////////////////// END OF : bidirectional graph SCC (deprecated and broken, dont use until proven otherwise)/////////////////////
///////////////////          (see start of this warning for notes)                            /////////////////////

// experimental: drop tips shared by multiple unitigs
void hamt_asgarc_ugTreatMultiLeaf(asg_t *sg, ma_ug_t *ug, int threshold_l){  // threshold_l is in bp, not number of reads or overlaps
    // FUNC
    //     remove arcs that have multiple arcs of the same direction on one side
    //     (aka try to resolve the web-like structure not much seen in other species)
    int verbose = 1;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv, v, w;
    asg_arc_t *av;
    int nb_del = 0, nb_del_node = 0, flag;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        flag = 0;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0)<2 || hamt_asgarc_util_countPre(auxsg, vu, 0, 0)>0){  // not the multi leaf
            continue;
        }
        if (ug->u.a[vu>>1].len>threshold_l){  // tip too long
            continue;
        }
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            wu = av[i].v;
            if (hamt_asgarc_util_isTip(auxsg, wu, 0, 0)){continue;}  // dont drop arc for tips
            hamt_ug_arc_del(sg, ug, vu, wu, 1);

            nb_del++;
            flag = 1;
        }
        if (flag){
            nb_del_node++;
            if (verbose){
                fprintf(stderr, "[debug::%s] trimmed utg%.6d\n", __func__, (int)(vu>>1)+1);
            }
        }
    }
    if (nb_del){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d links (%d nodes treated).\n", __func__, nb_del, nb_del_node);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}



// experimental: drop bifurcations that lead to different paths
//               and after dropping (dont call asg_cleanup!), 
//               rescue circles (since SCC is broken, can't efficiently/reliably avoid cutting circles)
void hamt_asgarc_ugTreatBifurcation(asg_t *sg, ma_ug_t *ug, int threshold_nread, int threshold_bp){
    int verbose = 1;
    double startTime = Get_T();
    asg_t *auxsg = ug->g;

    uint8_t *color1 = (uint8_t*)calloc(auxsg->n_seq*2, 1);  // buffer for branching 1
    uint8_t *color2 = (uint8_t*)calloc(auxsg->n_seq*2, 1);  // buffer for branching 2
    queue32_t q1, q2;
    queue32_init(&q1);
    queue32_init(&q2);
    uint32_t v0, v1, v2;  // unitig ID
    uint32_t va, vb, vc;  // vertex ID
    int nb_cut = 0, nb_rescue = 0, nb_skip = 0;
    int has_overlap = 0;

    // for rescuing
    uint8_t *color_sg = (uint8_t*)calloc(sg->n_seq*2, 1);
    int *subg_labels = (int*)calloc(sg->n_seq*2, sizeof(int));

    // cut bifurcations
    for (v0=0; v0<auxsg->n_seq*2; v0++){
        has_overlap = 0;
        if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0)!=2){
            continue;
        }
        int idx_v = 0;
        uint32_t tmp[2];
        asg_arc_t *tmpav = asg_arc_a(auxsg, v0);
        uint32_t tmpnv = asg_arc_n(auxsg, v0);
        for (uint32_t i=0; i<tmpnv; i++){
            if (!tmpav[i].del){
                tmp[idx_v] = tmpav[i].v;
                idx_v++;
            }
        }
        assert(idx_v==2);
        v1 = tmp[0];
        v2 = tmp[1];
        
        if (hamt_asgarc_util_countPre(auxsg, v1, 0, 0)!=1 || hamt_asgarc_util_countPre(auxsg, v2, 0, 0)!=1){  // require the bifucation to be simple
            continue;
        }

        // early termination: simple cases
        if (hamt_asgarc_util_checkSimpleBubble(auxsg, v0, NULL)){
            // v0 is the start of a simple bubble
            continue;
        }
        if (hamt_asgarc_util_countNoneTipSuc(auxsg, v1)<1 || hamt_asgarc_util_countNoneTipSuc(auxsg, v2)<1){
            // at least one side terminated in a tip
            continue;
        }

        ////// traversal //////
        // (not the most efficient way but this looks cleaner)
        // fprintf(stderr, "traversing utg%.6d dir %d.\n", (int)(v0>>1)+1, (int)(v0&1)); fflush(stderr);
        queue32_reset(&q1); queue32_reset(&q2);
        memset(color1, 0, auxsg->n_seq*2);  memset(color2, 0, auxsg->n_seq*2);
        if (hamt_ugasgarc_util_BFS_mark(ug, v1, color1, &q1, threshold_nread, threshold_bp)){  // early termination, v1 leads to a not-that-long path
            continue;
        }
        if (hamt_ugasgarc_util_BFS_mark(ug, v2, color2, &q2, threshold_nread, threshold_bp)){
            continue;
        }
        // (if we reach here, v1 and v2 goes to different & long paths, consider cutting the bifurcation)
        uint32_t san_i;
        for (uint32_t i1=0; i1<auxsg->n_seq*2; i1++){
            if (color1[i1] && color2[i1]){
                has_overlap = 1;
                san_i = i1;
            }
            if (has_overlap){
                break;
            }
        }
        if (has_overlap){  // the bifurcation doesn't lead to 2 diverging paths, dont cut
            // fprintf(stderr, "  overlap at %d.\n", (int)(san_i>>1)+1);
            continue;
        }

        // cut with circle sanchecks (the more careful way)
        // if (hamt_asgarc_detect_circleDFS2(auxsg, v0, v1)==0){
        //     hamt_ug_arc_del(sg, ug, v0, v1, 1);    
        //     nb_cut++;
        //     if (verbose){
        //         fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v0>>1)+1, (int)(v1>>1)+1);
        //     }
        // }else{
        //     nb_skip++;
        // }
        // if (hamt_asgarc_detect_circleDFS2(auxsg, v0, v0)==0){
        //     hamt_ug_arc_del(sg, ug, v0, v2, 1);
        //     nb_cut++;
        //     if (verbose){
        //         fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v0>>1)+1, (int)(v1>>1)+1);
        //     }
        // }else{
        //     nb_skip++;
        // }

        // just cut (the more aggressive way)
        hamt_ug_arc_del(sg, ug, v0, v1, 1);
        hamt_ug_arc_del(sg, ug, v0, v2, 1);
        fprintf(stderr, "the aggressive bifur cut: at utg%.6d, utg%.6d, utg%.6d\n", (int)(v0>>1)+1, (int)(v2>>1)+1, (int)(v2>>1)+1);
        nb_skip++;
        // // (exp: rescue at every cut)
        // // (note: nope dont do it, doesn't work at all.)
        // queue32_reset(&q1);  // repurpose
        // memset(color_sg, 0, sg->n_seq*2*sizeof(uint8_t));
        // memset(subg_labels, 0, sg->n_seq*2*sizeof(int));
        // int max_label = hamt_asgarc_util_BFS_markSubgraph(sg, subg_labels, color_sg, &q1); // note: this is on the string graph
        // int tmp_san;
        // for (int label=0; label<max_label; label++){
        //     nb_rescue += hamt_asgarc_rescueArcWithinCurrentSubgraph(sg, ug, v0);  // note: this is on the string graph
        //     nb_rescue += hamt_asgarc_rescueArcWithinCurrentSubgraph(sg, ug, v1);
        //     nb_rescue += hamt_asgarc_rescueArcWithinCurrentSubgraph(sg, ug, v2);
        // }    
    }

    // rescue: if there was an arc linking two tips of the same unconnected subgraph, recover it
    // TODO: maybe this is useless if we have the above circle detection
    //       it's also useless if a circle got more than 1 cut
    // uint8_t *color_sg = (uint8_t*)calloc(sg->n_seq*2, 1);
    // int *subg_labels = (int*)calloc(sg->n_seq*2, sizeof(int));
    // queue32_reset(&q1);  // repurpose
    // int max_label = hamt_asgarc_util_BFS_markSubgraph(sg, subg_labels, color_sg, &q1); // note: this is on the string graph
    // for (int label=0; label<max_label; label++){
    //     nb_rescue += hamt_asgarc_rescueArcWithinSubgraph(sg, NULL, subg_labels, label);  // note: this is on the string graph
    // }

    // clean up
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }

    queue32_destroy(&q1);
    queue32_destroy(&q2);
    free(color1);
    free(color2);
    free(subg_labels);
    free(color_sg);
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d arcs at bifurcation point, skipped %d; rescued %d arcs.\n", __func__, nb_cut, nb_skip, nb_rescue);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}

void hamt_ugasg_cut_shortTips(asg_t *sg, ma_ug_t *ug){
    // NOTE: only drop arcs, not deleting sequences
    int verbose = 1;

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    int pass = 1;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>30000){  // tip not short
            continue;
        }

        // specify the tip direction: no predecessor, has sucessor
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0)==0 || hamt_asgarc_util_countPre(auxsg, vu, 0, 0)!=0){
            continue;
        }

        // sancheck: don't drop tip if any of its sucessors has no target other than the tip
        pass = 1;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            wu = av[i].v;
            if (hamt_asgarc_util_countNoneTipSuc(auxsg, wu^1)==0){
                if (verbose){
                    fprintf(stderr, "[debug::%s] spared tip: utg%.6d \n", __func__, (int)(vu>>1)+1);
                }
                pass = 0;
                break;
            }
        }
        if (!pass){
            continue;
        }

        // cut
        for (uint32_t i=0; i<nv ;i++){
            if (av[i].del){continue;}
            wu = av[i].v;
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut tip: utg%.6d \n", __func__, (int)(vu>>1)+1);
            }
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d unitig tigs.\n", __func__, nb_cut);
    }
}

void hamt_asgarc_markBridges_DFS(asg_t *sg, uint32_t v, uint32_t v_parent, int *DFStime, uint8_t *visited, int *tin, int *low){
    // initended to be called (only) by hamt_asgarc_markBridges
    uint32_t w;
    asg_arc_t *av_fwd, *av_bwd, *av, *av_tmp;
    uint32_t nv_fwd, nv_bwd;
    visited[v>>1] = 1;
    tin[v>>1] = *DFStime;
    low[v>>1] = tin[v>>1];
    *DFStime = *DFStime+1;
    
    // check child nodes
    av_fwd = asg_arc_a(sg, v);
    av_bwd = asg_arc_a(sg, v^1);
    nv_fwd = asg_arc_n(sg, v);
    nv_bwd = asg_arc_n(sg, v^1);
    for (uint32_t i_=0, i=0; i_<(nv_fwd+nv_bwd); i_++){  // forward direction
        // fprintf(stderr, "at utg%.6d, time %d, i_=%d (n=%d)\n", (int)(v>>1)+1, *DFStime, (int)i_, (int)(nv_fwd+nv_bwd));
        if (i_<nv_fwd){
            i = i_;
            av = av_fwd;
        }else{
            i = i_-nv_fwd;
            av = av_bwd;
        }
        if (av[i].del){continue;}
        w = av[i].v;
        if (((v>>1)!=(v_parent>>1)) && (w>>1)==(v_parent>>1)){  // note: v==v_parent is the case for root node
            continue;
        }
        if (visited[w>>1]){
            low[v>>1] = low[v>>1]<tin[w>>1]? low[v>>1] : tin[w>>1];
        }else{
            hamt_asgarc_markBridges_DFS(sg, w, v, DFStime, visited, tin, low);
            low[v>>1] = low[v>>1]<low[w>>1]? low[v>>1] : low[w>>1];
            if (low[w>>1]>tin[v>>1]){  // arc is bridge
                av[i].is_bridge = 1;
                // mark the complementary arc as bridge too
                sg->arc[asg_arc_get_complementaryID(sg, &av[i])].is_bridge = 1;
            }
        }
    }
}

void hamt_ug_covCutByBridges(asg_t *sg, ma_ug_t *ug)
{
    // replacement of the bidirected SCC approach (which did not work as expected)
    // NOTE
    //     the graph is treated as if undirected, aka both directions will be touched when looking for child nodes
    int verbose = 1;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t v, w, nv;
    asg_arc_t *av;
    int cov_v, cov_w, cov_max, cov_min, nb_cut=0;
    
    // reset (just in case)
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        auxsg->arc[i].is_bridge = 0;
    }

    // mark bridge arcs
    uint8_t *visited = (uint8_t*)calloc(auxsg->n_seq, 1);
    int *tin = (int*)calloc(auxsg->n_seq, sizeof(int));  memset(tin, -1, auxsg->n_seq*sizeof(int));
    int *low = (int*)calloc(auxsg->n_seq, sizeof(int));  memset(tin, -1, auxsg->n_seq*sizeof(int));
    int DFStime = 0;
    for (v=0; v<auxsg->n_seq; v++){
        if (!visited[v>>1]){
            hamt_asgarc_markBridges_DFS(auxsg, v, v, &DFStime, visited, tin, low);
        }
    }
    if (verbose){
        int n_bridges = 0;
        for (uint32_t i=0; i<auxsg->n_arc; i++){
            if (auxsg->arc[i].is_bridge){
                n_bridges++;
                fprintf(stderr, "[debug::%s] bridge beween utg%.6d and utg%.6d\n", __func__, (int)(auxsg->arc[i].ul>>33)+1, (int)(auxsg->arc[i].v>>1)+1);
            }
        }
        fprintf(stderr, "[M::%s] %d out of %d arcs are bridge edges.\n", __func__, n_bridges, (int)auxsg->n_arc);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    // check arcs
    // (note: since the asg is treated as if undirected, edges in bubbles are not bridges)
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del){
            continue;
        }
        if (!auxsg->arc[i].is_bridge){
            continue;
        }
        
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        fprintf(stderr, "checking at utg%.6d vs utg%.6d\n", (int)(v>>1)+1, (int)(w>>1)+1);

        // check unitig coverage
        cov_v = ug->utg_coverage[v>>1];
        cov_w = ug->utg_coverage[w>>1];
        if (cov_v<cov_w){
            cov_min = cov_v;
            cov_max = cov_w;
        }else{
            cov_min = cov_w;
            cov_max = cov_v;
        }
        if (!hamt_check_covDiff(cov_min, cov_max)){
            fprintf(stderr, "failed cov diff: utg%.6d utg%.6d\n", (int)(v>>1)+1, (int)(w>>1)+1);
            continue;
        }

        // detect circles
        if (hamt_asgarc_detect_circleDFS(auxsg, v, w, 0)){  // note: v and w is ug with direction
            if (hamt_asgarc_util_isTip(auxsg, v, 0, 0) && (hamt_asgarc_util_countSuc(auxsg, v, 0, 0)==1)){
                // don't remove tips from circles
                // (however, allow cutting if the tip is still connected to somewhere else)
                continue;
            }

            // cut
            hamt_ug_arc_del(sg, ug, v, w, 1);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut: utg%.6d - utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                fprintf(stderr, "           (cov max: %d, cov min:%d)\n", cov_max, cov_min);
                fprintf(stderr, "             %d, %d; %d, %d\n", (int)hamt_asgarc_util_countPre(auxsg, v, 0, 0), 
                                                                 (int)hamt_asgarc_util_countSuc(auxsg, w, 0, 0), 
                                                                 (int)hamt_asgarc_util_countPre(auxsg, v, 0, 0), 
                                                                 (int)hamt_asgarc_util_countSuc(auxsg, w, 0, 0));
            }
        }else{
            fprintf(stderr, "   (failed circle test)\n");   
        }

        
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d links.\n", __func__, nb_cut);
    }
    free(tin);
    free(visited);
    free(low);
}

//////////////////////////////////////////////////////////
//                   contig gen                         //
//////////////////////////////////////////////////////////
void hamt_ug_pop_simpleSoapBubble(asg_t *sg, ma_ug_t *ug){
    // modified `drop_semi_circle`
    // pop at the stem, not middle of the soap bubble

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu, nvu;
    asg_arc_t *avu;

    for (vu=0; vu<auxsg->n_seq*2;vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0)==0 | hamt_asgarc_util_countPre(auxsg, vu, 0, 0)==0){
            continue;
        }
        // seek to the left


        // seek to the right


        // is a soap bubble, remove arcs joining the stem

        


    }

}



//////////////////////////////////////////////////////////
//                   interface                          //
//////////////////////////////////////////////////////////


void hamt_circle_cleaning(asg_t *sg, ma_ug_t *ug)
{
    // hamt_asgarc_ugCovCutSCC(sg, ug, coverage_cut, sources, ruIndex);  
    hamt_ug_covCutByBridges(sg, ug);
    hamt_asgarc_ugCovCutDFSCircle(sg, ug);
}

void hamt_clean_shared_seq(asg_t *sg, ma_ug_t *ug){
    hamt_asgarc_ugTreatMultiLeaf(sg, ug, 50000);
    hamt_asgarc_ugTreatBifurcation(sg, ug, -1, 10000000);
}

void hamt_ug_pop_bubble(asg_t *sg,ma_ug_t *ug)
{
    // remove simple bubbles at unitig graph level
    // NOTE: only drop arcs, not deleting sequences
    int verbose = 1;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu[2], uu, nv;
    asg_arc_t *av;
    // uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));  // for utg coverage
    int cov[2], idx, nb_cut=0;

    // pre clean
    hamt_ugasg_cut_shortTips(sg, ug);

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0)!=2){
            continue;
        }
        if (!hamt_asgarc_util_checkSimpleBubble(auxsg, vu, &uu)){  // check if the current utg is the start of a simple bubble
            continue;
        }
        av = asg_arc_a(auxsg, vu);
        for (int i=0, i_=0; i<asg_arc_n(auxsg, vu); i++){
            if (av[i].del) {continue;}
            wu[i_++] = av[i].v;
        }
        // check unitig length
        if (ug->u.a[wu[0]>>1].len>50000 || ug->u.a[wu[1]>>1].len>50000){
            continue;
        }

        // check coverage, retain the edge with higher coverage    
        // cov[0] = get_ug_coverage(&ug->u.a[wu[0]>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
        // cov[1] = get_ug_coverage(&ug->u.a[wu[1]>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
        cov[0] = ug->utg_coverage[wu[0]>>1];
        cov[1] = ug->utg_coverage[wu[1]>>1];
        idx = cov[0]>cov[1]? 1:0;
        // drop arc
        if (verbose){
            fprintf(stderr, "[debug::%s]> utg%.6d\n", __func__, (int)(vu>>1)+1);
            fprintf(stderr, "             vu is utg%.6d, wu is %.6d, uu is %.6d\n", (int)(vu>>1)+1, (int)(wu[idx]>>1)+1, (int)(uu>>1)+1);
        }
        hamt_ug_arc_del(sg, ug, vu, wu[idx], 1);
        hamt_ug_arc_del(sg, ug, wu[idx], uu, 1);
        nb_cut++;
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d bubbles on unitig graph.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    // free(primary_flag);
}

void hamt_ug_pop_miscbubble(asg_t *sg, ma_ug_t *ug){
    // drop alternative edges, even if it's not a simple bubble
    // affected topo should include:
    //   - multi-edge "bubbles"
    //   - not-strictly-a-bubble bubbles (i.e. start/end node has other linkage)
    int verbose = 1;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu[20], uu[20];
    uint32_t nv;
    int idx;
    int color[20], flag, nb_cut = 0;
    asg_arc_t *av;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>50000){
            continue;
        }

        idx = 0;
        memset(color, 0, sizeof(int)*20);
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0)<3){
            continue;
        }
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){
                continue;
            }
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0)!=1 || hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0)!=1){
                continue;
            }
            if (ug->u.a[av[i].v>>1].len>50000){
                continue;
            }
            wu[idx] = av[i].v;
            hamt_asgarc_util_get_the_one_target(auxsg, av[i].v, &uu[idx], 0, 0);
            idx++;
            if (idx==20){  // buffer sancheck
                fprintf(stderr, "[W::%s] more than 20 targets?? \n", __func__);
                break;
            }
        }
        if (idx==0){
            continue;
        }

        // check if any edge shares uu with another edge; preserve the 1st one and discard all others
        for (int i1=0; i1<idx; i1++){
            if (color[i1]){continue;}  // this wu has been cut
            if (hamt_asgarc_util_isTip(auxsg, uu[i1], 0, 0)){continue;}  // ignore misc bubble that ends at a tig
            flag = 0;
            for (int i2=i1; i2<idx; i2++){
                if (color[i2]){continue;}
                if (uu[i2]!=uu[i1]){continue;}
                color[i2] = 1;
                if (flag){  // cut
                    hamt_ug_arc_del(sg, ug, vu, wu[i2], 1);
                    hamt_ug_arc_del(sg, ug, wu[i2], uu[i2], 1);
                    if (verbose){
                        fprintf(stderr, "[debug::%s] cut, u=utg%.6d, w=utg%.6d, u=utg%.6d\n", __func__,
                                          (int)(vu>>1)+1, (int)(wu[i2]>>1)+1, (int)(uu[i2]>>1)+1);
                    }
                    nb_cut++;
                }else{
                    flag = 1;
                }
            }
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d spots\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}

int hamt_ug_pop_simpleInvertBubble(asg_t *sg, ma_ug_t *ug){
    int verbose = 1;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu;
    int nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>50000){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleInvertBubble(sg, vu, &wu, &uu)){
            if (ug->u.a[wu>>1].len>50000){
                continue;
            }
            if (ug->u.a[uu>>1].len>50000){
                continue;
            }
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            hamt_ug_arc_del(sg, ug, wu, uu, 1);
            nb_cut++;

            if (verbose){
                fprintf(stderr, "[debug::%s] start=utg%.6d, cut=utg%.6d, end=utg%.6d \n", __func__, 
                                (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
            }
        }
    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d locations\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut;
}

int hamt_sg_pop_simpleInvertBubble(asg_t *sg){
    //  used by the hamt preclean
    int verbose = 1;

    uint32_t v, w, u;
    int nb_cut = 0;

    for (v=0; v<sg->n_seq*2; v++){
        if (hamt_asgarc_util_checkSimpleInvertBubble(sg, v, &w, &u)){
            asg_arc_del(sg, v, w, 1);
            asg_arc_del(sg, w^1, v^1, 1);
            asg_arc_del(sg, w, u, 1);
            asg_arc_del(sg, u^1, w^1, 1);
            nb_cut++;
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
    }
    return nb_cut;
}

void hamt_ug_prectgTopoClean(asg_t *sg){
    // note: will redo the unitig graph
    // assume graph is clean
    // pop all simple bubbles, apparent bubbles, pentagon bubbles (popping doesn't depend on coverage)
    double startTime = Get_T();

    ma_ug_t *ug = ma_ug_gen(sg);
    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu, nv;
    asg_arc_t *av;
    int nb_pop = 0, round = 0;

    for (round=0; round<5; round++){
        nb_pop = 0;
        // 2-edge or multi-edge simple bubbles
        for (vu=0; vu<auxsg->n_seq*2; vu++){
            if (!hamt_asgarc_util_checkSimpleBubble_multiEddge(auxsg, vu, &uu)){
                continue;
            }
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            int san = 0;
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del){continue;}
                if (hamt_asgarc_util_isTip(auxsg, av[i].v, 0, 0)){continue;}
                wu = av[i].v;
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                hamt_ug_arc_del(sg, ug, wu, uu, 1);    
                san++;
                // fprintf(stderr, "ugpopdebug\tutg%.6d, utg%.6d, utg%.6d\n", (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
            }
            // keep the last edge
            hamt_ug_arc_del(sg, ug, vu, wu, 0);
            hamt_ug_arc_del(sg, ug, wu, uu, 0);
            if (san>1){
                // TODO: bug or not, this is probably caused by allowing tips and something was off
                nb_pop++;
            }
        }
        hamtdebug_add_fake_utg_coverage(ug);
        hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, 990+round);
        free(ug->utg_coverage);
        fprintf(stderr, "[M::%s] popped %d simple bubbles (round %d).\n", __func__, nb_pop, round);

        // inverted simple bubbles
        // nb_pop += hamt_ug_pop_simpleInvertBubble(sg, ug);

        if (nb_pop==0){
            ma_ug_destroy(ug);
            if (VERBOSE){
                fprintf(stderr, "[M::%s] (early termination at round %d)\n", __func__, round);
            }
            break;
        }else{
            asg_cleanup(sg);
            ma_ug_destroy(ug); ug = ma_ug_gen(sg); auxsg = ug->g;
            if (VERBOSE){
                fprintf(stderr, "[M::%s] popped total %d locations (round %d).\n", __func__, nb_pop, round);
                fprintf(stderr, "[T::%s] took %0.2fs.\n\n", __func__, Get_T()-startTime);

            }
        }
    }

}

void hamt_ug_prectg_rescueShortCircuit(asg_t *sg, ma_ug_t *ug){
    // exp
    // addresses observations in the sheep datasets, where 2 apparent circles
    //   are joined by one relatively short unitig (not necessarily 1-read)
    //   if the ends of the circle actually overlap with each other, then ditch the short unitig and 
    //   recover that arc (which was dropped during transitive reduction)
    // FUNC
    //    recover some arc
    // NOTE
    //    will regenerate the unitig graph, maybe several times
    //    hopefully this doesn't mess up with other stuffs?
    



}