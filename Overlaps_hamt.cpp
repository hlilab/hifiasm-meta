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

#define HAMT_PRIMARY_LABEL 0
#define HAMT_ALTER_LABEL 1

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
    // simplified `ma_ug_print2`
    // NOTE
    //    Assume coverage has been collected
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

void hamtdebug_output_unitig_graph_ug(ma_ug_t *ug, char *base_fn, const char *suffix, int cleanID){
    // write gfa without ma_ug_gen and sequence
    // solely meant for graph cleaning debug
    char* gfa_name = (char*)malloc(strlen(base_fn)+100);
    sprintf(gfa_name, "%s.G%s_%d.utg.gfa", base_fn, suffix, cleanID);
    FILE *output_file = fopen(gfa_name, "w");
    hamt_ug_print_simple(ug, &R_INF, output_file);
    fclose(output_file);
    free(gfa_name);
}

//////////////////////////////////////////////////////////////////////
//                 non-debug file output routines                   //
//////////////////////////////////////////////////////////////////////
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
    if ((q->n+2)>=q->m){
        q->a = (uint32_t*)realloc(q->a, sizeof(uint32_t)*(q->m + (q->m>>1)));
        assert(q->a);
        q->prev_m = q->m;
        q->m = q->m + (q->m>>1);
        if (q->i_head>q->i_tail){
            q->earlywrap = 1;
        }
    }
    q->a[q->i_tail] = d;
    q->n++;
    if (q->i_tail==(q->m-1)){
        q->i_tail = 0;
    }else{
        q->i_tail++;
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

int queue32_get_size(queue32_t *q){
    return q->n;
}

int queue32_is_in_queue(queue32_t *q, uint32_t d){
    // check if a given value is in the current queue
    // RETURN
    //     0 if no
    //     1 if yes
    if (q->n==0){return 0;}  // queue is empty
    int i = q->i_head;
    while (i!=q->i_tail){
        if (q->a[i]==d){return 1;}
        if (!q->earlywrap){
            if (i==(q->m-1)){
                i = 0;
            }else{
                i++;
            }
        }else{
            if (i==(q->prev_m-1)){
                i = 0;
            }else{
                i++;
            }
        }
    }
    if (q->a[i]==d){return 1;}
    return 0;
}


typedef struct {
    uint64_t *a;
    int n, m;
}vecu64_t;
void vecu64_init(vecu64_t *v){
    v->n = 0;
    v->m = 16;
    v->a = (uint64_t*)malloc(16*sizeof(uint64_t));
}
void vecu64_destroy(vecu64_t *v){
    free(v->a);
}
void vecu64_push(vecu64_t *v, uint64_t d){
    if (v->n==v->m){
        v->m = v->m + (v->m>>1);
        v->a = (uint64_t*)realloc(v->a, sizeof(uint64_t)*v->m);
    }
    v->a[v->n++] = d;
}
int vecu64_is_in_vec(vecu64_t *v, uint64_t d){
    for (int i=0; i<v->n; i++){
        if (v->a[i]==d){
            return 1;
        }
    }
    return 0;
}
int vecu64_get_size(vecu64_t *v){
    return v->n;
}

typedef struct {
    uint32_t *a;
    int n, m;
}vecu32_t;
void vecu32_init(vecu32_t *v){
    v->n = 0;
    v->m = 16;
    v->a = (uint32_t*)malloc(16*sizeof(uint32_t));
}
void vecu32_destroy(vecu32_t *v){
    free(v->a);
}
void vecu32_push(vecu32_t *v, uint32_t d){
    if (v->n==v->m){
        v->m = v->m + (v->m>>1);
        v->a = (uint32_t*)realloc(v->a, sizeof(uint32_t)*v->m);
    }
    v->a[v->n++] = d;
}
void vecu32_reset(vecu32_t *v){
    v->n = 0;
}
int vecu32_is_in_vec(vecu32_t *v, uint32_t d){
    for (int i=0; i<v->n; i++){
        if (v->a[i]==d){
            return 1;
        }
    }
    return 0;
}
int vecu32_get_size(vecu32_t *v){
    return v->n;
}


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

int stack32_is_in_stack(stacku32_t *stack, uint32_t d){
    // FUNC
    //     check if a given value is in the buffer
    // RETURN
    //     1 if yes
    //     0 if no
    if (stack->n==0){return 0;}
    for (int i=0; i<stack->n; i++){
        if (stack->a[i]==d){
            return 1;
        }
    }
    return 0;
}





//////////////////////////////////////////////////////////////////////////
//                        routines for overlaps                         //
//////////////////////////////////////////////////////////////////////////

int hamt_ovlp_read_coverage_nbreads(ma_hit_t_alloc *sources, long long i_read, uint32_t qs, uint32_t qe){
    // FUNC
    //     Roughly estimate read ocverage by counting overlap lengths of the given read.
    if (qe<=qs){return 0;}
    float ret=0; 
    ma_hit_t *h;
    uint32_t s, e;
    for (int i=0; i<sources[i_read].length; i++){
        h = &sources[i_read].buffer[i];
        if (h->del) {continue;}
        s = (uint32_t)h->qns > qs? (uint32_t)h->qns : qs;
        e = h->qe < qe? h->qe : qe;
        if (e>s) {ret+=e-s;}
    }
    ret = ret/(qe-qs);
    return ret<=1? 1 : (int)ret;  // waterproof zero for easier use
}


//////////////////////////////////////////////////////////////////////
//                        helper routines                           //
//////////////////////////////////////////////////////////////////////

int does_ovlp_ever_exist(ma_hit_t_alloc *sources, uint32_t v0, uint32_t w0, ma_hit_t **h0, ma_hit_t **h0_rev){
    // NOTE
    //     sources is effectively R_INF.paf
    // FUNC
    //     check overlaps of the same haplotype
    // RETURN
    //     1 if yes
    //     0 if no
    int verbose = 0;

    int nb_targets;
    int found = 0;
    // v0 to w0
    nb_targets = sources[v0>>1].length;
    for (int i=0; i<nb_targets; i++){
        if (sources[v0>>1].buffer[i].tn==(w0>>1)){
            if (verbose){fprintf(stderr, "[debug::%s] found (direction 1)\n", __func__);}
            found++;
            *h0 = &sources[v0>>1].buffer[i];
            break;
        }
    }
    // w0^1 to v0^1
    nb_targets = sources[w0>>1].length;
    for (int i=0; i<nb_targets; i++){
        if (sources[w0>>1].buffer[i].tn==(v0>>1)){
            if (verbose){fprintf(stderr, "[debug::%s] found (direction 2)\n", __func__);}
            found++;
            *h0_rev = &sources[w0>>1].buffer[i];
            break;
        }
    }
    if (found==2){return 1;}
    if (verbose){fprintf(stderr, "[debug::%s] failed\n", __func__);}
    return 0;
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
    fprintf(stderr, "- droparc: utg%.6d - utg%.6d\n", (int)(vu>>1)+1, (int)(wu>>1)+1);
}

static inline void hamt_ug_utg_softdel(asg_t *sg, ma_ug_t *ug, uint32_t vu, int label){
    // label unitig vu and its corresponding reads
    uint32_t v;
    for (int i=0; i<ug->u.a[vu>>1].n; i++){
        v = (uint32_t) (ug->u.a[vu>>1].a[i]>>33);
        sg->seq[v].c = label;
    }
    ug->u.a[vu>>1].c = label;
    ug->g->seq_vis[vu>>1] = label;  // note seq_vis is actually n_seq*2 long
    if (ug->u.a[vu>>1].len>0){  // debug
        fprintf(stderr, "- softdel: utg%.6d (length %d), to label %d\n", (int)(vu>>1)+1, (int)ug->u.a[vu>>1].len, label);
    }
}

int hamt_ug_arc_del_between(asg_t *sg, ma_ug_t *ug, uint32_t start, uint32_t end, int del, int label){
    // treat every arc between vu and wu
    // NOTE
    //     set label to -1 to ignore sequence label
    uint8_t *color = (uint8_t*)calloc(ug->g->n_seq*2, 1);
    uint32_t vu, nv;
    asg_arc_t *av;
    asg_t *auxsg = ug->g;
    int nb_cut = 0;
    int yes_terminate = 0;

    queue32_t q;
    queue32_init(&q);
    queue32_enqueue(&q, start);

    while (queue32_dequeue(&q, &vu)){
        if (color[vu]==2) {continue;}
        color[vu] = 2;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].v==end){
                yes_terminate = 1;
                // fprintf(stderr, "debug: terminating... finish remaining stuff\n");
                // break;
            }
            if (color[av[i].v]!=0) {continue;}
            color[av[i].v] = 1;
            if (!yes_terminate){
                queue32_enqueue(&q, av[i].v);
            }
            hamt_ug_arc_del(sg, ug, vu, av[i].v, del);
            nb_cut++;
            // fprintf(stderr, "debug: cut_betwee - utg%.6d , utg%.6d\n", (int)(vu>>1)+1, (int)(av[i].v>>1)+1);
        }
        // if (yes_terminate){
        //     break;
        // }
    }
    free(color);
    queue32_destroy(&q);
    return nb_cut;
}

int hamt_ug_arc_flip_between(asg_t *sg, ma_ug_t *ug, uint32_t start, uint32_t end, int base_label, int new_label){
    // mark every vertex between start and end with new_label
    // NOTE
    //     set label to -1 to ignore sequence label
    uint8_t *color = (uint8_t*)calloc(ug->g->n_seq*2, 1);
    uint32_t vu, nv;
    asg_arc_t *av;
    asg_t *auxsg = ug->g;
    int nb_cut = 0;
    int yes_terminate = 0;

    queue32_t q;
    queue32_init(&q);
    queue32_enqueue(&q, start);

    while (queue32_dequeue(&q, &vu)){
        if (color[vu]==2) {continue;}
        color[vu] = 2;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            // if (ug->u.a[vu>>1].c!=base_label) {continue;}
            if (av[i].v==end){
                yes_terminate = 1;
                // fprintf(stderr, "debug: terminating... finish remaining stuff\n");
                // break;
            }
            if (color[av[i].v]!=0) {continue;}
            color[av[i].v] = 1;
            if (!yes_terminate){
                queue32_enqueue(&q, av[i].v);
            }
            // mark
            if (av[i].v!=end){
                hamt_ug_utg_softdel(sg, ug, av[i].v, new_label);
            }
            nb_cut++;
        }
        // if (yes_terminate){
        //     break;
        // }
    }
    free(color);
    queue32_destroy(&q);
    return nb_cut;
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

void hamt_ug_cleanup_arc_by_labels(asg_t *sg, ma_ug_t *ug){
    // NOTE
    //     ma_ug_gen_primary only looks at the label when seeding unitig,
    //     so for soft cut we still need to drop some arcs.
    asg_t *auxsg = ug->g;
    uint32_t vu, wu;
    int cnt = 0;
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del) {continue;}
        vu = (uint32_t) (auxsg->arc[i].ul>>32);
        wu = auxsg->arc[i].v;
        if (auxsg->seq_vis[vu>>1]!=auxsg->seq_vis[wu>>1]){
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            cnt++;
        }
    }
    if (cnt){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[debug::%s] cleaned up %d arcs\n", __func__, cnt);
    }
}

void hamt_ug_init_seq_vis(ma_ug_t *ug, uint8_t flag){
    if (!ug->g->seq_vis){
        ug->g->seq_vis = (uint8_t*)calloc(ug->g->n_seq*2, 1);
    }else{
        memset(ug->g->seq_vis, 0, ug->g->n_seq*2);
    }
}
ma_ug_t *hamt_ug_gen(asg_t *sg,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex, int flag){
    ma_ug_t *ug;
    if (flag<0){
        ug = ma_ug_gen(sg);
    }else{
        ug = ma_ug_gen_primary(sg, flag);
    }
    hamt_ug_init_seq_vis(ug, flag);
    hamt_collect_utg_coverage(sg, ug, coverage_cut, sources, ruIndex);
    return ug;
}
void hamt_ug_destroy(ma_ug_t *ug){
    hamt_destroy_utg_coverage(ug);
    ma_ug_destroy(ug);
}
void hamt_ug_regen(asg_t *sg, ma_ug_t **ug,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex, int flag){
    hamt_ug_cleanup_arc_by_labels(sg, *ug);
    hamt_ug_destroy(*ug);
    *ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, flag);
}
void hamt_asg_reset_seq_label(asg_t *sg, uint8_t flag){
    int primary=0, alter=0;
    for (uint32_t i=0; i<sg->n_seq; i++){
        if (sg->seq[i].c==HAMT_PRIMARY_LABEL){primary++;}
        else if (sg->seq[i].c==HAMT_ALTER_LABEL){alter++;}
        sg->seq[i].c = flag;
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] reset seq labels. (Stat before: primary %d, alter %d)\n", __func__, primary, alter);
    }
}


void hamt_asg_reset_seq_vis(asg_t *sg, uint8_t flag){
    // NOTE
    //    hamt recycles seq_vis allocation for auxsg (i.e. ug->g)
    //    for safety concern, it might be desireble for the sg
    for (uint32_t i=0; i<sg->n_seq; i++){
        sg->seq_vis[i] = flag;  // seq_vis is actually of length n_seq*2
    }
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

int hamt_check_diploid(ma_ug_t *ug, uint32_t vu1, uint32_t vu2, float ratio0,
                       ma_hit_t_alloc *reverse_sources)
{
    // FUNC
    //     (similar to check_if_diploid)
    //     Check if vu1 and vu2 appears to be haplotigs.
    // NOTE
    //     Doesn't do topo check for vu1 and vu2. 
    //     The caller is responsible for checking topo context.
    // RETURN
    //     -1 if no hit at all
    //     1 if yes
    //     0 if no
    asg_t *auxsg = ug->g;
    uint32_t vu_short, vu_long;
    uint32_t qn, qn2, tn;
    int nb_targets;
    int found, cnt_total=0, cnt_hit=0;
    float ratio = ratio0<0? 0.3 : ratio0;
    
    if (ug->u.a[vu1>>1].len < ug->u.a[vu2>>1].len){
        vu_short = vu1;
        vu_long = vu2;
    }else{
        vu_short = vu2;
        vu_long = vu1;
    }

    for (int i=0; i<ug->u.a[vu_short>>1].n; i++){  // iterate over reads in the shorter unitig
        qn = ug->u.a[vu_short>>1].a[i]>>33;

        // check if there's any read in the other unitig that targets the read 
        found = 0;
        for (int j=0; j<ug->u.a[vu_long>>1].n; j++){  // iterate over reads in the other unitig
            qn2 = ug->u.a[vu_long>>1].a[j]>>33;
            nb_targets = reverse_sources[qn2].length;
            for (int i_target=0; i_target<nb_targets; i_target++) {  // check all existed inter-haplotype targets of this read
                if (reverse_sources[qn2].buffer[i_target].tn == qn){
                    found = 1;
                    break;
                }
            }
            if (found){break;}
        }

        // update stats
        cnt_total++;
        if (found){
            cnt_hit++;
        }
    }
    if (cnt_hit==0){return -1;}
    if ((float)cnt_hit/cnt_total > ratio ){return 1;}
    return 0;
}

/////////////////////////
//         asg         //
/////////////////////////

static inline int hamt_asgarc_util_countPre(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
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
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        nv++;
    }
    return nv;
}
static inline int hamt_asgarc_util_countSuc(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
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
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        nv++;
    }
    return nv;
}

int hamt_asgarc_util_isTip(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
    // FUNC
    //     check if v is a tip (or orphan), in either direction
    // RETURN
    //     0 if no
    //     1 if yes
    if (hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc, base_label)==0 || 
        hamt_asgarc_util_countPre(g, v, include_del_seq, include_del_arc, base_label)==0){
        return 1;
    }
    return 0;
}

int hamt_asgarc_util_isDanglingTip(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
    // FUNC
    //     check if v is a tip (or orphan), in either direction
    // RETURN
    //     0 if no
    //     1 if yes
    if (hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc, base_label)==1 || 
        hamt_asgarc_util_countPre(g, v, include_del_seq, include_del_arc, base_label)==0){
        return 1;
    }else if (hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc, base_label)==0 || 
              hamt_asgarc_util_countPre(g, v, include_del_seq, include_del_arc, base_label)==1){
        return 1;
    }
    return 0;
}


int hamt_asgutil_detect_strict_single_long_path(asg_t *sg, uint32_t begNode, uint32_t *endNode, 
                                                int threshold_nodes, uint64_t threshold_length, int base_label){
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
        nv = hamt_asgarc_util_countSuc(sg, v, 0, 0, base_label);
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
        if (hamt_asgarc_util_countPre(sg, v, 0, 0, base_label)!=1) return 0;  // backward bifur 
    }
    return 1;
}


int hamt_asgarc_util_get_the_one_target(asg_t *g, uint32_t v, uint32_t *w, int include_del_seq, int include_del_arc, int base_label){
    // v as exactly one target (sucessor), given the `include_del_seq` and `include_del_arc` (NOT tested in this function)
    // get that vertex's index in av (NOT vertex ID)
    int nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    for (i=0; i<(int)nv; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        if (w){
            *w = av[i].v;
        }
        return i;
    }
    // fprintf(stderr, "[E::%s] unexpected.\n", __func__); exit(1);
    return -1;
}
int hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(asg_t *g, uint32_t v, uint32_t *w, 
                                                          int include_del_seq, int include_del_arc, int base_label){
    // v as exactly one target (sucessor), given the `include_del_seq` and `include_del_arc` (NOT tested in this function)
    // get that vertex's index in av (NOT vertex ID)
    int nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    for (i=0; i<(int)nv; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        if (hamt_asgarc_util_countSuc(g, av[i].v, 0, 0, base_label)==0 && hamt_asgarc_util_countPre(g, av[i].v, 0, 0, base_label)==1){continue;}
        if (w){
            *w = av[i].v;
        }
        return i;
    }
    // fprintf(stderr, "[E::%s] unexpected.\n", __func__); exit(1);
    return -1;
}

int hamt_asgarc_util_walk_furthest_strict_single_path(asg_t *g, uint32_t v0, uint32_t *vn, int *nb_nodes, int limit, int base_label){
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
        if (hamt_asgarc_util_countSuc(g, v, 0, 0, base_label)!=1){
            break;
        }
        hamt_asgarc_util_get_the_one_target(g, v, &v, 0, 0, base_label);
        if (hamt_asgarc_util_countPre(g, v, 0, 0, base_label)>0){
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

int hamt_asgarc_util_countSinglePath(asg_t *g, uint32_t v0, int include_del_seq, int include_del_arc, int base_label, 
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
        nb_fwd = hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc, base_label);
        if (nb_fwd==0) {return -2;}  // reaches the end of this tip
        if (nb_fwd>1) {break;}  // the end, or more than 1 targets
        wid = hamt_asgarc_util_get_the_one_target(g, v, &w, include_del_seq, include_del_arc, base_label);
        av = asg_arc_a(g, v);
        if (wid==-1){
            fprintf(stderr, "waht?\n");fflush(stderr);
            exit(1);
        }
        if ((w>>1) == (v0>>1)){
            return -1;  // loop
        }

        nb_bwd =  hamt_asgarc_util_countPre(g, w, include_del_seq, include_del_arc, base_label);
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

int hamt_asgarc_util_countNoneTipSuc(asg_t *g, uint32_t v, int base_label){
    // NOTE
    //    actually, no dangling tip
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
        if (hamt_asgarc_util_countSuc(g, w, 0, 0, base_label)==0){
            continue;
        }
        cnt++;
    }
    return cnt;
}
int hamt_asgarc_util_countNoneDanglingTipSuc(asg_t *g, uint32_t v, int base_label){
    // NOTE
    //    actually, no dangling tip
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
        if (hamt_asgarc_util_countSuc(g, w, 0, 0, base_label)==0 && hamt_asgarc_util_countPre(g, w, 0, 0, base_label)==1){
            continue;
        }
        cnt++;
    }
    return cnt;
}


int hamt_asgarc_util_checkSimpleBubble(asg_t *g, uint32_t v0, uint32_t *dest, int base_label){
    // NOTE
    //     tolerates 1-vertex tip
    // FUNC
    //     check if v0 is the start of a simple bubble (i.e. a 1-vertex bubble)
    // RETURN
    //     0 if no
    //     1 if yes
    //     -1 if circular
    int verbose = 0;
    // if (hamt_asgarc_util_countSuc(g, v0, 0, 0)!=2){
    if (hamt_asgarc_util_countNoneDanglingTipSuc(g, v0, base_label)!=2){
        if (verbose){fprintf(stderr, "[checkBubble] v-suc\n");}
        return 0;
    }
    uint32_t w[2], nv, u, u1, u2, nv_tmp;  // x is the next-next vertex (if exists)
    asg_arc_t *av = asg_arc_a(g, v0), *av_tmp;
    int idx = 0;
    int idx_target;
    
    for (uint32_t i=0; i<asg_arc_n(g, v0); i++){
        if (av[i].del) {continue;}
        if (g->seq[av[i].v>>1].del) {continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        // if (hamt_asgarc_util_isDanglingTip(g, av[i].v, 0, 0)){continue;}
        if (hamt_asgarc_util_countSuc(g, av[i].v, 0, 0, base_label)==0 && hamt_asgarc_util_countPre(g, av[i].v, 0, 0, base_label)==1){continue;}

        // // important: don't allow a circle to be a simple bubble!
        // //            (check the target of the persumed bubble edge)
        // nv_tmp = asg_arc_n(g, av[i].v);
        // av_tmp = asg_arc_a(g, av[i].v);
        // for (uint32_t i_tmp=0; i_tmp<nv_tmp; i_tmp++){
        //     if ((av_tmp[i].v>>1)==(v0>>1)){
        //         return -1;
        //     }
        // }

        w[idx] = av[i].v;
        idx++;
    }
    // assert(idx==2);
    if (idx!=2){
        fprintf(stderr, "assertion failed, idx is %d\n", idx);
        exit(1);
    }
    if (verbose){
        fprintf(stderr, "  (w1 utg%.6d w2 utg%.6d)\n", (int)(w[0]>>1)+1, (int)(w[1]>>1)+1);
    }

    // check if the next-next vertex is the end of a bubble
    // if (hamt_asgarc_util_countSuc(g, w[0], 0, 0)!=1 || hamt_asgarc_util_countSuc(g, w[1], 0, 0)!=1){  // need to end in only one vertex
    if (hamt_asgarc_util_countNoneTipSuc(g, w[0], base_label)!=1 || hamt_asgarc_util_countNoneTipSuc(g, w[1], base_label)!=1){  // need to end in only one vertex
        if (verbose){fprintf(stderr, "[checkBubble] w-suc\n");}
        return 0;
    }
    if (hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(g, w[0], &u1, 0, 0, base_label)==-1){fprintf(stderr, "[%s] ??????", __func__);return 0;}
    if (hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(g, w[1], &u2, 0, 0, base_label)==-1){fprintf(stderr, "[%s] ??????", __func__);return 0;}
    if (u1!=u2){  // and it's the same vertex
        if (verbose){fprintf(stderr, "[checkBubble] u\n");}
        return 0;
    }
    // if (hamt_asgarc_util_countPre(g, u1, 0, 0)!=2){  // and the "end of the simple bubble" must not have other incoming sources
    if (hamt_asgarc_util_countNoneTipSuc(g, u1^1, base_label)!=2){  // and the "end of the simple bubble" must not have other incoming sources
        if (verbose){fprintf(stderr, "[checkBubble] u-pre (u is utg%.6d , pre: %d)\n", (int)(u1>>1)+1, hamt_asgarc_util_countPre(g, u1, 0, 0, base_label));}
        return 0;
    }
    // check backward links of the bubble edges
    // if (hamt_asgarc_util_countPre(g, w[0], 0, 0)!=1 || hamt_asgarc_util_countPre(g, w[1], 0, 0)!=1){
    if (hamt_asgarc_util_countNoneTipSuc(g, w[0]^1, base_label)!=1 || hamt_asgarc_util_countNoneTipSuc(g, w[1]^1, base_label)!=1){
        if (verbose){fprintf(stderr, "[checkBubble] w-suc\n");}
        return 0;
    }
    if (dest){
        *dest = u1;
    }
    if ((u1>>1)!=(v0>>1)){
        return 1;
    }else{
        return -1;
    }
}

int hamt_asgarc_util_checkSimpleBubble_edge(asg_t *g, uint32_t v0, int base_label){
    // check if v0 is an edge of a simple bubble (i.e. 1-vertex bubble)
    //     (doesn't allow tip, but could've; modify if needed)
    // RETURN
    //     0 if no
    //     1 if yes

    uint32_t s, e, e2, sib;  // start, end, sibling
    uint32_t nv;
    asg_arc_t *av;

    // check simple edge
    if (hamt_asgarc_util_countPre(g, v0, 0, 0, base_label)!=1 || hamt_asgarc_util_countSuc(g, v0, 0, 0, base_label)!=1){
        return 0;
    }
    // check both ends of the persumed bubble
    if (hamt_asgarc_util_get_the_one_target(g, v0, &e, 0, 0, base_label)<0){
        fprintf(stderr, "[ERROR::%s] e\n", __func__);fflush(stderr);
        exit(1);
    }
    if (hamt_asgarc_util_get_the_one_target(g, v0^1, &s, 0, 0, base_label)<0){
        fprintf(stderr, "[ERROR::%s] s\n", __func__);fflush(stderr);
        exit(1);
    }
    s = s^1;
    if (hamt_asgarc_util_countPre(g, e, 0, 0, base_label)!=2 || hamt_asgarc_util_countSuc(g, s, 0, 0, base_label)!=2){
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
    if (hamt_asgarc_util_countPre(g, sib, 0, 0, base_label)!=1 || hamt_asgarc_util_countSuc(g, sib, 0, 0, base_label)!=1){  // compile note: is safe
        return 0;
    }
    if (hamt_asgarc_util_get_the_one_target(g, sib, &e2, 0, 0, base_label)<0){
        fprintf(stderr, "[ERROR::%s] e2\n", __func__);fflush(stderr);
        exit(1);
    }
    if (e2==e){
        return 1;
    }
    return 0;
}

int hamt_asgarc_util_checkSimpleBubble_multiEddge(asg_t *sg, uint32_t v0, uint32_t *u0, int base_label){
    // check if v0 is the start of a multi-edge simple bubble (including 2-edge simple bubble)
    // RETURN
    //      0 if no
    //      1 if yes (also right the end vertex to *u0)
    uint32_t w, u_tmp, u, nv, nw;
    int idx = 0;
    asg_arc_t *av, *aw;
    if (hamt_asgarc_util_countNoneTipSuc(sg, v0, base_label)<2){
        return 0;
    }
    nv = asg_arc_n(sg, v0);
    av = asg_arc_a(sg, v0);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        if (hamt_asgarc_util_isTip(sg, av[i].v, 0, 0, base_label)){continue;}
        w = av[i].v;
        if (hamt_asgarc_util_countPre(sg, w, 0, 0, base_label)!=1 || hamt_asgarc_util_countNoneTipSuc(sg, w, base_label)!=1){
            return 0;
        }
        // check each edge's target
        aw = asg_arc_a(sg, w);
        for (uint32_t i2=0; i2<asg_arc_n(sg, w); i2++){  // find the non-tip target
            if (aw[i2].del){continue;}
            if (hamt_asgarc_util_isTip(sg, aw[i2].v, 0, 0, base_label)){continue;}
            u_tmp = aw[i2].v;
            if (idx==0){
                u = u_tmp;
            }else{
                if (u_tmp!=u){  // compile note: this is safe
                    return 0;
                }
            }
            idx++;
            break;
        }
    }
    if (idx==0){
        return 0;
    }
    if (idx!=hamt_asgarc_util_countNoneTipSuc(sg, u^1, base_label)){
        return 0;
    }
    if ((u>>1)==(v0>>1)){
        // TODO: forgot why i need this, but i do
        return 0;
    }
    *u0 = u;
    return 1;
}

int hamt_asgarc_util_checkSimpleInvertBubble(asg_t *sg, uint32_t v0, uint32_t *w0, uint32_t *u0, int base_label){
    // FUNC
    /*
            ----w1--
          /         \
    .....v0         |
          \        /
           -----u1-----...
          we want to drop v0->w1 and w1->u1^1
    */
   //     only treat the most simple case, any branching and it'll be ignored
   // RETURN
   //     1 if v0 is the startpoint of such structure
   //     0 otherwise
    int verbose = 0;

    uint32_t uw[2], x;
    int nuw[4];
    uint32_t nv, idx=0;
    asg_arc_t *av;
    if (hamt_asgarc_util_countSuc(sg, v0, 0, 0, base_label)!=2){
        return 0;
    }
    nv = asg_arc_n(sg, v0);
    av = asg_arc_a(sg, v0);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        uw[idx] = av[i].v;
        idx++;
    }
    assert(idx==2);
    
    // check topo
    nuw[0] = hamt_asgarc_util_countSuc(sg, uw[0], 0, 0, base_label);
    nuw[1] = hamt_asgarc_util_countSuc(sg, uw[1], 0, 0, base_label);
    nuw[2] = hamt_asgarc_util_countPre(sg, uw[0], 0, 0, base_label);
    nuw[3] = hamt_asgarc_util_countPre(sg, uw[1], 0, 0, base_label);
    if (verbose){fprintf(stderr, "[debug::%s] at ID %.6d\n", __func__, (int)(v0>>1)+1);}
    if (nuw[0]==0 || nuw[1]==0){  // no target
        if (verbose){fprintf(stderr, "[debug::%s]    no target\n", __func__);}
        return 0;
    }
    if (nuw[2]>1 || nuw[3]>1){  // branching in backward direction
        if (verbose){fprintf(stderr, "[debug::%s]    backward branching\n", __func__);}
        return 0;
    }
    // one and only one vertex shall have two targets
    int idx_rev;
    if (nuw[0]==1 && nuw[1]==2){
        idx_rev = 0;  // uw[0] is the w1 in the topo comment
    }else if (nuw[0]==2 && nuw[1]==1){
        idx_rev = 1;  // uw[1] is the w1
    }else{
        if (verbose){fprintf(stderr, "[debug::%s]    failed 1-and-2\n", __func__);}
        return 0;
    }
    // check topo: the inversion
    if (hamt_asgarc_util_get_the_one_target(sg, uw[idx_rev], &x, 0, 0, base_label)<0){
        fprintf(stderr, "[W::%s] failed to get target when supposed to.\n", __func__);fflush(stderr);
        exit(1);
    }
    if (x!=(uw[idx_rev^1]^1)){
        if (verbose){fprintf(stderr, "[debug::%s]    not revserse link\n", __func__);}
        return 0;
    }else{
        *w0 = uw[idx_rev];
        *u0 = uw[idx_rev^1]^1;
        if (verbose){fprintf(stderr, "[debug::%s]    CUT\n", __func__);}
        return 1;
    }

}

int hamt_asgarc_util_checkSimplePentaBubble(asg_t *g, uint32_t v0, uint32_t *buffer, int base_label){
    // FUNC
    //    pops the following structure:
    /*       v2---v5.
           /   \    \ 
       ---v1    v4  v7---
           \     \ /
            v3---v6
        (or the mid link be v3->v4->v5 instead of v2->v4->v6. it's the same.)
    */
   //     by removing the 2 links (4 arcs) in path v2-v4-v6 
   // NOTE
   //    buffer is at least length 3.
   //    doesn't tolerate tips.
   //    also assumes only arc can be marked as del, not sequences (like most other funcs in hamt graph cleaning)
   // RETURN
   //    1 if found
   //    0 if none
    int verbose = 0;
    uint32_t w1[2], w2[2], w3[2], w_tmp;
    uint32_t v2, v3, v4, v5, v6, v7, v7_;
    uint32_t nv;
    asg_arc_t *av;
    int idx;

    if (hamt_asgarc_util_countSuc(g, v0, 0, 0, base_label)!=2){
        if (verbose){fprintf(stderr, "    v0 failed\n");}
        return 0;
    }
    // step 1st vertex
    nv = asg_arc_n(g, v0);
    av = asg_arc_a(g, v0);
    idx = 0;
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label) {continue;}
        if (hamt_asgarc_util_countPre(g, av[i].v, 0, 0, base_label)!=1){
            if (verbose){fprintf(stderr, "    w1 failed\n");}
            return 0;
        }
        w1[idx] = av[i].v;
        idx++;
    }
    if (idx!=2){
        fprintf(stderr, "[W::%s] abnormal idx (is %d), continue anyway\n", __func__, idx);
        return 0;
    }
    // check step 1

    if (hamt_asgarc_util_countSuc(g, w1[0], 0, 0, base_label)==1 && hamt_asgarc_util_countSuc(g, w1[1], 0, 0, base_label)==2){
        v2 = w1[1];
        v3 = w1[0];
    }else if (hamt_asgarc_util_countSuc(g, w1[0], 0, 0, base_label)==2 && hamt_asgarc_util_countSuc(g, w1[1], 0, 0, base_label)==1){
        v2 = w1[0];
        v3 = w1[1];
    }else{
        if (verbose){fprintf(stderr, "    w1 suc failed\n");}
        return 0;
    }

    // step the single edge (v3->v6)
    if (hamt_asgarc_util_get_the_one_target(g, v3, &v6, 0, 0, base_label)<0){
        fprintf(stderr, "[W::%s] failed at v3->v6, continue anyway\n", __func__);
        return 0;
    }
    if (hamt_asgarc_util_countPre(g, v6, 0, 0, base_label)!=2 || hamt_asgarc_util_countSuc(g, v6, 0, 0, base_label)!=1){
        if (verbose){fprintf(stderr, "    v6 topo failed\n");}
        return 0;
    }

    // check v2->v4->v6
    nv = asg_arc_n(g, v2);
    av = asg_arc_a(g, v2);
    idx = 0;
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label) {continue;}
        if (hamt_asgarc_util_countSuc(g, av[i].v, 0, 0, base_label)!=1 || hamt_asgarc_util_countPre(g, av[i].v, 0, 0, base_label)!=1){
            if (verbose){fprintf(stderr, "    w2 failed\n");}
            return 0;
        } 
        w2[idx] = av[i].v;
        if (hamt_asgarc_util_get_the_one_target(g, av[i].v, &w3[idx], 0, 0, base_label)<0){
            fprintf(stderr, "[W::%s] error at getting w2's target, continue anyway\n", __func__);
            return 0;
        }
        idx++;
    }
    if (idx!=2){
        fprintf(stderr, "[W::%s] 2nd abnormal idx (is %d), continue anyway\n", __func__, idx);
        return 0;
    }
    if (w3[0]==v6){
        v4 = w2[0];
        v5 = w2[1];
        v7 = w3[1];
    }else if (w3[1]==v6){
        v4 = w2[1];
        v5 = w2[0];
        v7 = w3[0];
    }else{
        if (verbose){fprintf(stderr, "    w3 failed\n");}
        return 0;
    }

    // check sink vertex (v5 and v6's shared target)
    if (hamt_asgarc_util_get_the_one_target(g, v6, &v7_, 0, 0, base_label)<0){
        fprintf(stderr, "[W::%s] failed to get v6's target, continue anyway\n", __func__);
        return 0;
    }
    if (v7!=v7_){
        if (verbose){fprintf(stderr, "    sink failed\n");}
        return 0 ; 
    }
    
    buffer[0] = v2;
    buffer[1] = v4;
    buffer[2] = v6;
    return 1;
}

int hamt_ug_util_popSimpleBiBubbleChain(asg_t *sg, ma_ug_t *ug, uint32_t v0, 
                                        int base_label){
    // FUNC
    //    pops the following structure 
    //    (this can be chained indefinitely long, as long as it fall into a v6-like sink eventually, or ends up in tips):
    /*       v2---v4
           /   \   \ 
       ---v1   |  v6---
           \   \  /
            v3-v5
    */
   //     by removing all v2-v5-like links
   // NOTE
   //    always prompt the caller to remove the shorter side
   //    doesn't tolerate tips.
   //    also assumes only arc can be marked as del, not sequences (like most other funcs in hamt graph cleaning)
   // RETURN
   //    number of cut performed
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t w1[2], w2[2], w3[2], w_tmp;
    uint32_t v2, v3, v4, v5, v6, v7, v7_;
    uint32_t nv;
    asg_arc_t *av;
    int idx, color[2], n_suc, n_pre, linked;
    int nb_cut = 0;

    uint32_t prv_va=0, prv_vb=0;

    // check the handle
    if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0, base_label)!=2){
        if (verbose){fprintf(stderr, "    v0 failed\n");}
        return 0;
    }
    nv = asg_arc_n(auxsg, v0);
    av = asg_arc_a(auxsg, v0);
    idx = 0; color[0] = 0; color[1] = 0;
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)>1){return 0;}
        n_suc = hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label);
        if (n_suc>2 || n_suc==0){return 0;}
        // the vertex has either 1 or 2 targets
        if (av[i].v==v0){return 0;}  // paranoia
        w1[n_suc-1] = av[i].v;
        color[n_suc-1] = 1;
        idx++;
    }
    assert(idx==2);
    // require that one and only one of {v2, v3} to have 1 target, and the other one must have 2
    if (color[0]==0 || color[1]==0){return 0;}

    // early termination case: only one bi-bubble
    uint32_t u1;
    linked = 0;
    int passed = 1;
    if (hamt_asgarc_util_get_the_one_target(auxsg, w1[0], &u1, 0, 0, base_label)<0){
        fprintf(stderr, "ERROR %s u1\n", __func__);
        exit(1);
    }
    if (hamt_asgarc_util_countSuc(auxsg, u1, 0, 0, base_label)<=1){
        nv = asg_arc_n(auxsg, w1[1]);
        av = asg_arc_a(auxsg, w1[1]);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (av[i].v==u1){
                linked = 1;
            }
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label)>1){passed = 0;}
        }
        if (linked && passed){
            // cut 
            hamt_ug_arc_del(sg, ug, w1[1], u1, 1);
            if (verbose){fprintf(stderr, "    DROPPED utg%.6d - utg%.6d\n", (int)(prv_va>>1)+1, (int)(prv_vb>>1)+1);}
            nb_cut++;
            return nb_cut;
        }
    }


    // walk the persumed bubble chain
    if (verbose){fprintf(stderr, "    initial check ok, enter walk\n");}
    // while (1){
    int san;
    for (san=0; san<20; san++){
        // get the target in the single edge, sancheck topo
        if (hamt_asgarc_util_get_the_one_target(auxsg, w1[0], &w_tmp, 0, 0, base_label)<0){
            fprintf(stderr, "error1 %s\n", __func__);
            exit(1);
        }
        if (w_tmp==v0){return nb_cut;}  // paranoia
        if (hamt_asgarc_util_countPre(auxsg, w_tmp, 0, 0, base_label)!=2){
            return nb_cut;
        }
        if (hamt_asgarc_util_countSuc(auxsg, w_tmp, 0, 0, base_label)>2){
            return nb_cut;
        }
        if (verbose){fprintf(stderr, "    checkpoint 1\n");}

        // get the two target on the other side, sancheck topo
        nv = asg_arc_n(auxsg, w1[1]);
        av = asg_arc_a(auxsg, w1[1]);
        idx = 0; color[0] = 0; color[1] = 0; linked = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)>2){return nb_cut;}
            n_suc = hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label);
            n_pre = hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label);
            if (n_suc>2){return nb_cut;}
            if (w_tmp==v0){return nb_cut;}  // paranoia

            w2[n_pre-1] = av[i].v;
            color[n_pre-1] = 1;
            idx++;
        }
        assert(idx==2);
        if (color[0]==0 || color[1]==0){
            return nb_cut;
        }
        if (verbose){fprintf(stderr, "    checkpoint 2\n");}

        // check if there's a cross arc
        if (w_tmp!=w2[1]){return nb_cut;}
        if (verbose){fprintf(stderr, "    checkpoint 3\n");}

        // check w2's targets
        if (hamt_asgarc_util_countSuc(auxsg, w2[0], 0, 0, base_label)==1 && hamt_asgarc_util_countSuc(auxsg, w2[1], 0, 0, base_label)==2){
            w3[0] = w2[0];
            w3[1] = w2[1];
        }else if (hamt_asgarc_util_countSuc(auxsg, w2[0], 0, 0, base_label)==2 && hamt_asgarc_util_countSuc(auxsg, w2[1], 0, 0, base_label)==1){
            w3[0] = w2[1];
            w3[1] = w2[0];
        }else{
            // chain terminated
            // check the end is a sink or two tips
            //    if yes, cut previously logged arc (if exist)
            //    otherwise, do nothing
            // TODO/BUG
            //    this part itself will miss the last shouldve-been-cut arc of a chain
            //    however this function is meant to be called multiple times, so practically ok right now. (Nov 5 2020)
            //    Come back later and fix it.
            int suc1 = hamt_asgarc_util_countSuc(auxsg, w2[0], 0, 0, base_label);
            int suc2 = hamt_asgarc_util_countSuc(auxsg, w2[1], 0, 0, base_label);
            uint32_t tmp_v1, tmp_v2;
            if (suc1==1 && suc2==1){
                if (hamt_asgarc_util_get_the_one_target(auxsg, w2[0], &tmp_v1, 0, 0, base_label)<0){
                    fprintf(stderr, "ERROR %s, tmp_v1\n", __func__);
                    exit(1);
                }
                if (hamt_asgarc_util_get_the_one_target(auxsg, w2[1], &tmp_v2, 0, 0, base_label)<0){
                    fprintf(stderr, "ERROR %s, tmp_v2\n", __func__);
                    exit(1);
                }
                if (hamt_asgarc_util_countPre(auxsg, tmp_v1, 0, 0, base_label)<=2 && hamt_asgarc_util_countPre(auxsg, tmp_v2, 0, 0, base_label)<=2){
                    if (!(prv_va==0 && prv_vb==0)){
                        hamt_ug_arc_del(sg, ug, prv_va, prv_vb, 1);
                        if (verbose){fprintf(stderr, "    DROPPED utg%.6d - utg%.6d\n", (int)(prv_va>>1)+1, (int)(prv_vb>>1)+1);}
                        nb_cut++;
                    }
                }
            }
            return nb_cut;
        }
        if (verbose){fprintf(stderr, "    checkpoint 4\n");}

        // current topo checks passed, cut previous one and log the current one
        // (cut)
        if (prv_va==0 && prv_vb==0){  // first round, init
            prv_va = w1[1];
            prv_vb = w2[1];  // aka w_tmp
        }else{
            hamt_ug_arc_del(sg, ug, prv_va, prv_vb, 1);
            if (verbose){fprintf(stderr, "    DROPPED utg%.6d - utg%.6d\n", (int)(prv_va>>1)+1, (int)(prv_vb>>1)+1);}
            nb_cut++;
        }
        // (update)
        if (verbose){
            fprintf(stderr, "    checkpoint 5, prv_va utg%.6d prv_vb utg%.6d\n", (int)(prv_va>>1)+1, (int)(prv_vb>>1)+1);
        }
        prv_va = w1[1];
        prv_vb = w2[1];
        w1[0] = w3[0];
        w1[1] = w3[1];

    }
    if (san==20){
        fprintf(stderr, "HEY PROBABLY INFINITE LOOP\n");fflush(stderr);
        exit(1);
    }
    return nb_cut;
}


int hamt_asgarc_util_countSinglePath_allowSimpleBubble(asg_t *g, uint32_t v0, int max_arcs, int base_label){
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
        nb_fwd = hamt_asgarc_util_countSuc(g, v, 0, 0, base_label);
        if (nb_fwd==0) {  // reached the end of this single path
            return -1;
        } else if (nb_fwd==1){
            hamt_asgarc_util_get_the_one_target(g, v, &w, 0, 0, base_label);
            nb_bwd = hamt_asgarc_util_countPre(g, w, 0, 0, base_label);
            if (nb_bwd!=1){
                return 1;  // the next target has other incoming sources, terminate
            }else{  // step
                v = w;
                dist++;
            }
        } else if (nb_fwd==2){
            int ret = hamt_asgarc_util_checkSimpleBubble(g, v, &w, base_label);
            if (ret>0){
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

int hamt_ugasg_util_findendSinglePath_allowSimpleBubble(ma_ug_t *ug, uint32_t v0, uint32_t *end0, int *is_circle, int *length_bp, int base_label){
    // FUNC
    //    given v0, walk in both directions and find the end of the single path
    //    and store them in *end0
    //    (directions are: v0->...->end0)
    // RETURN
    //    length (in terms of vertices; 1 bubble is counted as 1)
    // NOTE
    //    tolerate 1-vertex tip (if g is sg, then 1-read tip; if g is auxsg, it's a 1-utg tip)
    int verbose = 0;
    int cnt = 0;
    asg_t *g = ug->g;
    uint32_t v, w, u, nv;
    asg_arc_t *av;
    v = v0;
    *end0 = v;
    uint8_t *color = (uint8_t*)calloc(g->n_seq, 1);
    color[v>>1] = 1;

    *is_circle = 0;

    if (verbose){
        fprintf(stderr, "[debug::%s] at v0 utg%.6d\n", __func__, (int)(v0>>1)+1);
    }
    while (1){
        if (verbose){ fprintf(stderr, "[debug::%s]  v is utg%.6d\n", __func__, (int)(v>>1)+1); }
        if ((hamt_asgarc_util_checkSimpleBubble(g, v^1, NULL, base_label))<0){  // circular structure found when identifying simple bubble
            *is_circle = 1;
            if (verbose){ fprintf(stderr, "[debug::%s]     leave (detected circle)\n", __func__); }
            break;
        }
        if (hamt_asgarc_util_countNoneDanglingTipSuc(g, v^1, base_label)!=1 && (hamt_asgarc_util_checkSimpleBubble(g, v^1, NULL, base_label))<=0){ 
            // the vertex has multiple incomes, and it's not a end-of-simple-bubble edge
            if (verbose){ fprintf(stderr, "[debug::%s]     leave (backward branching)\n", __func__); }
            break;
        }
        if (hamt_asgarc_util_countNoneDanglingTipSuc(g, v, base_label)==0){  // end of walk
            cnt++;
            *end0 = v;
            if (verbose){ fprintf(stderr, "[debug::%s]     leave (end of walk: utg%.6d)\n", __func__, (int)(v>>1)+1); }
            break;
        }

        if (hamt_asgarc_util_countNoneDanglingTipSuc(g, v, base_label)>1){ // step through a bubble
            if (verbose){ fprintf(stderr, "[debug::%s]     check if start of simple bubble\n", __func__); }
            if (hamt_asgarc_util_checkSimpleBubble(g, v, &u, base_label)>0){
                // sancheck: circling
                // if end of the bubble loops back to a previously seen vertex, don't step this bubble
                if (color[u>>1]){  
                    *is_circle = 1;
                    if (verbose){ fprintf(stderr, "[debug::%s]     break (v0 in the circle)\n", __func__); }
                    break;
                }

                // update
                cnt+=2;
                *end0 = v;
                v = u;
                color[u>>1] = 1;
                if (verbose){ fprintf(stderr, "[debug::%s]     stepped bubble\n", __func__); }
            }else{  // non-simple-bubble branching
                if (verbose){ fprintf(stderr, "[debug::%s]     leave (forward non-bubble branching)\n", __func__); }
                break;
            }
        }else{  // step to the next vertex
            if (hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(g, v, &w, 0, 0, base_label)<0){
                fprintf(stderr, "ERROR AT %s\n", __func__);
                exit(1);
            }
            if (color[w>>1]){  
                *is_circle = 1;
                if (verbose){ fprintf(stderr, "[debug::%s]     break (v0 in the circle)\n", __func__); }
                break;
            }

            // TODO/NOTE: length of simple bubble's edge isn't counted
            //            because when writting this function, i only needed a rough estimation
            //            aka the path to be not too short, but the threshold itself is already arbitrary.
            //            If in the future we need the exact length, modify it.
            *length_bp += ug->u.a[v>>1].len;

            cnt+=1;
            *end0 = v;
            v = w;
            color[w>>1] = 1;
            if (verbose){ fprintf(stderr, "[debug::%s]     step to next vertex (utg%.6d)\n", __func__, (int)(w>>1)+1); }
        }
    }
    free(color);
    return cnt;
}

int hamt_ugasgarc_util_BFS_mark(ma_ug_t *ug, uint32_t v0, uint8_t *color, queue32_t *q,
                                int threshold_nread, int threshold_bp, int base_label){
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
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
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


int hamt_asgarc_util_BFS_markSubgraph_core(asg_t *sg, int *subg_labels, uint8_t *color0, queue32_t *q0, int base_label){
    // helper func, do a BFS within the given range(s) and mark touched vertices in `color`
    // (if both thresholds are given, require exceed them both to terminate search)
    // (set threshold to -1 to disable)
    // the caller is responsible for init color and q
    // RETURN
    //      max subgraph label
    // TODO
    //     subg_labels and color don't need to be n_seq*2 long. 
    int verbose = 0;

    uint32_t v0, v, nv, w;
    asg_arc_t *av;
    int cnt_nread, cnt_bp;
    int tmp;
    int label = 0;

    queue32_t *q = 0;
    queue32_t q_alloc;
    uint8_t *color;
    if (!q0){
        queue32_init(&q_alloc);
        q = &q_alloc;
    }else{
        q = q0;
    }
    if (!color0){
        color = (uint8_t*)calloc(sg->n_seq*2, 1);
    }else{
        color = color0;
    }

    for (v0=0; v0<sg->n_seq*2; v0++){
        if (color[v0]){
            continue;
        }
        if (base_label>=0 && sg->seq_vis[v0>>1]!=base_label){continue;}
        tmp = 0;
        queue32_reset(q);
        queue32_enqueue(q, v0);
        color[v0] = 1;
        color[v0^1] = 1;
        tmp++;
        subg_labels[v0] = label;
        subg_labels[v0^1] = label;
        while(queue32_dequeue(q, &v)){
            if (color[v]==2){
                continue;
            }
            if (base_label>=0 && sg->seq_vis[v>>1]!=base_label){continue;}
            // 1st direction
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            for (uint32_t i=0; i<nv; i++){
                // if (av[i].del){continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(q, w);
                    color[w] = 1;
                    color[w^1] = 1;
                    tmp++;
                    subg_labels[w] = label;
                    subg_labels[w^1] = label;
                }
            }
            // 2nd direction
            nv = asg_arc_n(sg, v^1);
            av = asg_arc_a(sg, v^1);
            for (uint32_t i=0; i<nv; i++){
                // if (av[i].del){continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(q, w);
                    color[w] = 1;
                    color[w^1] = 1;
                    tmp++;
                    subg_labels[w] = label;
                    subg_labels[w^1] = label;
                }
            }
            color[v] = 2;
            color[v^1] = 2;
        }
        label++;
        if (verbose){
            fprintf(stderr, "[debug::%s] subgraph #%d, %d reads.\n", __func__, label, tmp);
        }
    }

    if (!q0){queue32_destroy(&q_alloc);}
    if (!color0){free(color);}
    return label;
}

void hamt_ug_util_BFS_markSubgraph(ma_ug_t *ug, int base_label){
    asg_t *auxsg = ug->g;
    uint8_t *color_sg = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    int *subg_labels = (int*)calloc(auxsg->n_seq*2, sizeof(int));
    queue32_t q;
    queue32_init(&q);
    int max_label = hamt_asgarc_util_BFS_markSubgraph_core(auxsg, subg_labels, color_sg, &q, base_label);
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){
        ug->u.a[i>>1].subg_label = subg_labels[i];
    }
    free(color_sg);
    free(subg_labels);
    queue32_destroy(&q);
}

int hamt_asgutil_calc_singlepath_cov(asg_t *g, uint32_t v, int base_label){
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
        if (hamt_asgarc_util_countSuc(g, w, 0, 0, base_label)!=1){
            break;
        }
        wid = hamt_asgarc_util_get_the_one_target(g, w, &w, 0, 0, base_label);
        // w = asg_arc_a(g, w)[wid].v;
        if (hamt_asgarc_util_countPre(g, w, 0, 0, base_label)!=1){
            break;
        }
    }
    cov = cov/n;
    return (int) (cov+0.499);
}

int hamt_asgutil_is_tigSafeToErode(asg_t *g, uint32_t v, int base_label){
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
    has_branching = hamt_asgarc_util_countSinglePath_allowSimpleBubble(g, v, 10, base_label);  // the number of vertices, not bp
    if (has_branching==1){
        return 1;  // this tip can be erode
    }else{
        return 0;  // don't touch it
    }

}

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

int hamt_sg_pop_simpleInvertBubble(asg_t *sg){
    int verbose = 0;

    uint32_t v, w, u;
    int nb_cut = 0;

    for (v=0; v<sg->n_seq*2; v++){
        if (hamt_asgarc_util_checkSimpleInvertBubble(sg, v, &w, &u, -1)){
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

void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length){
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


int hamt_asgarc_detect_circleDFS(asg_t *sg, uint32_t v0, uint32_t w0, int allow_v0, int base_label){ // vu0 has direction; v2vu has directions
    // !experimental notes!
    //     on tolerating inversions: 
    //       when scaning a node, if we found a white target whose other direction is gray
    //       push the edge^1, mark the white target as black; only push target that is white in both direction
    //     TODO
    //       handle checking is probably not always correct
    // FUNC
    //     DFS on sg (could be ug's auxsg), test if vu is the start of a circle
    //     (note "start of a circle" means v0 is the "stem"; v0 shall not be part of the circle)
    // RETURN
    //     return 1 if yes, 0 otherwise

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
            if (hamt_asgarc_util_countSuc(sg, v, 0, 0, base_label)==0){
                color[v] = 2;
                continue;
            }
            int is_exhausted = 1;
            int special_flag = 0;
            for (int i=0; i<(int)nv; i++){
                if (av[i].del){continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
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
                    // inversion heuristics
                    if (color[w^1]==1){
                        special_flag = 1;
                        color[w] = 2;  // attempt to not traverse back; TODO: might be buggy
                        continue;
                    }
                    // end of inversion heuristics

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
            if (special_flag){
                // current node found a white target that's gray in the other direction
                // therefore try the other direction for the current node too
                stacku32_push(&stack, v^1);
            }
        }else{
            break;
        }
    }
    free(color);
    stacku32_destroy(&stack);
    return ret;
}

int hamt_asgarc_detect_circleDFS2(asg_t *sg, uint32_t v0, uint32_t w0, uint32_t *v_exclude, int base_label){
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
                if (av[i].del) {continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label){continue;}
                if (color[av[i].v]==0){
                    is_exhausted = 0;
                    break;
                }
            }
            if (!is_exhausted){
                stacku32_push(&stack, v);
                for (uint32_t i=0; i<nv; i++){
                    if (av[i].del) {continue;}
                    if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label){continue;}
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
                if (av[i].del) {continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label){continue;}
                if (color[av[i].v]==0){
                    is_exhausted = 0;
                    break;
                }
            }
            if (!is_exhausted){
                stacku32_push(&stack, v);
                for (uint32_t i=0; i<nv; i++){
                    if (av[i].del) {continue;}
                    if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label){continue;}
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

void hamt_asgarc_ugCovCutDFSCircle(asg_t *sg, ma_ug_t *ug, int base_label)
{
    int verbose = 0;  // debug
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t v, w, nv;
    uint32_t cov[3];
    uint32_t covtmp;
    asg_arc_t *av;

    int nb_cut =0;

    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del){continue;}
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (base_label>=0 && (auxsg->seq_vis[v>>1]!=base_label || auxsg->seq_vis[w>>1]!=base_label) ) {continue;}
        if (hamt_asgarc_util_countPre(auxsg, w, 0, 0, base_label)<2){
            continue;
        }

        cov[2] = 0;
        av = asg_arc_a(auxsg, w^1);
        for (uint32_t i=0; i<asg_arc_n(auxsg, w^1); i++){  // collect the largest coverage of w's predecessor (except v which will be accounted by cov[0])
            if ((av[i].v>>1)==(v>>1)){  // ignore v
                continue;
            }
            if (av[i].del){continue;}
            if (base_label>=0 && (auxsg->seq_vis[av[i].v>>1]!=base_label)) {continue;}
            // covtmp = get_ug_coverage(&ug->u.a[av[i].v>>1], sg, coverage_cut, sources, ruIndex, primary_flag);
            covtmp = ug->utg_coverage[av[i].v>>1];
            if (covtmp>cov[2]){
                cov[2] = covtmp;
            }
        }
        cov[0] = ug->utg_coverage[v>>1];
        cov[1] = ug->utg_coverage[w>>1];
        radix_sort_ovhamt32(cov, cov+3);

        // check coverage
        if (!hamt_check_covDiff(cov[0], cov[2])){
            if (verbose>1) {fprintf(stderr, "[debug::%s] failed covDiff utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);}
            continue;
        }

        // detect circles and cut
        if (hamt_asgarc_detect_circleDFS(auxsg, v, w, 0, base_label)){  // note: v and w is ug with direction
            hamt_ug_arc_del(sg, ug, v, w, 1);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
            }
        }
    }
    
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d links.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
}


void hamt_asgarc_ugCovCutDFSCircle_aggressive(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //      when performing DFS-based circle detection,
    //      dont require the handle to be not included in the circle
    //      criteria
    //         1) there exists a circling path, without touching the handle
    //         2) the circling path involving the handle shall be much longer than the one in 1)
    //         3) more relaxed pre/suc count requirements and coverage diff requirements
    // NOTE
    //     experimental, only use it after all more carefule steps have been done
    int verbose = 0;  // debug
    double startTime = Get_T();
    int nb_cut = 0;

    asg_t *auxsg = ug->g;
    uint32_t v, w, nv, nw;
    asg_arc_t *av, *aw;
    uint32_t covs[20], cov_min, cov_max, cov_v, cov_tmp;  // TODO: tho it's probably safe enough

    for (v=0; v<auxsg->n_seq*2; v++){
        if (base_label>=0 && sg->seq_vis[v>>1]!=base_label) {continue;}
        // skip tips and simple bubble starts
        if (hamt_asgarc_util_countNoneTipSuc(auxsg, v, base_label)==0){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleBubble(auxsg, v, NULL, base_label)>0){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleBubble_edge(auxsg, v, base_label)>0){
            continue;
        }

        av = asg_arc_a(auxsg, v);
        nv = asg_arc_n(auxsg, v);
        for (int iw=0; iw<nv; iw++){
            if (av[iw].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[iw].v>>1]!=base_label) {continue;}
            w = av[iw].v;
            if (hamt_asgarc_util_isTip(auxsg, w, 0, 0, base_label)){
                continue;
            }
            if (hamt_asgarc_util_countPre(auxsg, w, 0, 0, base_label)<2){
                continue;
            }
            // check coverage diff
            nw = asg_arc_n(auxsg, w^1);
            aw = asg_arc_a(auxsg, w^1);
            int idx = 0;
            for (int i=0; i<nw; i++){
                if (aw[i].del){continue;}
                if (base_label>=0 && auxsg->seq_vis[aw[i].v>>1]!=base_label) {continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (hamt_asgarc_util_isTip(auxsg, aw[i].v, 0, 0, base_label)){continue;}  // don't use tip's coverage
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
                if (base_label>=0 && auxsg->seq_vis[aw[i].v>>1]!=base_label) {continue;}
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
                if (base_label>=0 && auxsg->seq_vis[aw[i].v>>1]!=base_label) {continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (hamt_asgarc_detect_circleDFS2(auxsg, w^1, aw[i].v, &v, base_label)>0){
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
                if (verbose>1){
                    fprintf(stderr, "[debug::%s] DIDNOT cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
            }
        }
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    }
    asg_cleanup(sg);
    asg_cleanup(auxsg);

}


int hamt_asgarc_rescueArcWithinSubgraph(asg_t *sg, ma_ug_t *ug, int *subg_labels, int label, int base_label){
    // (experimental) recover previous deleted arcs within the subgraph 
    // NOTE
    //     (this is not the function that checks arcs reduced by transitive redcution and recover them
    //       this function only reverts arc[i].del status for vertices of the same subgraph.)
    //     subg_labels: labeling of not connected subgraphs 
    //     label      : specify subgraph
    //     base_label : haplotype label (-1 to disable)
    uint32_t v, w, vu, wu;
    int cnt = 0;
    if (!ug){
        for (uint32_t i=0; i<sg->n_arc; i++){
            if (!sg->arc[i].del){
                continue;
            }
            v = sg->arc[i].ul>>32;
            w = sg->arc[i].v;
            if (base_label>=0 && (sg->seq_vis[v>>1]!=base_label || sg->seq_vis[w>>1]!=base_label)) {continue;}
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
            if (base_label>=0 && (auxsg->seq_vis[vu>>1]!=base_label || auxsg->seq_vis[wu>>1]!=base_label)) {continue;}
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

int hamt_asgarc_rescueArcWithinCurrentSubgraph(asg_t *sg, ma_ug_t *ug, uint32_t v0, int base_label){
    // dont need to label all subgraphs, just search based on the given seed
    int verbose = 0;
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
        if (base_label>=0 && (auxsg->seq_vis[v>>1]!=base_label || auxsg->seq_vis[w>>1]!=base_label)) {continue;}
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
void hamt_asgarc_ugTreatMultiLeaf(asg_t *sg, ma_ug_t *ug, int threshold_l,// threshold_l is in bp, not number of reads or overlaps
                                   int base_label, int alt_label, int is_hard_drop){  
    // FUNC
    //     remove arcs that have multiple arcs of the same direction on one side
    //     (aka try to resolve the web-like structure not much seen in other species)
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv, v, w;
    asg_arc_t *av;
    int nb_del = 0, nb_del_node = 0, flag;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        flag = 0;
        // multileaf tip direction: no predecessor, more than 2 sucessors
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<3 || hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)>0){
            continue;
        }
        // don't touch if it's a long tip
        if (ug->u.a[vu>>1].len>threshold_l){
            continue;
        }

        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            wu = av[i].v;
            // if (hamt_asgarc_util_isTip(auxsg, wu, 0, 0)){continue;}  // dont drop arc for tips
            if (is_hard_drop){
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
            }else{
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                // hamt_ug_utg_softdel(sg, ug, vu, alt_label);
            }

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
void hamt_asgarc_ugTreatBifurcation(asg_t *sg, ma_ug_t *ug, int threshold_nread, int threshold_bp, int base_label){
    int verbose = 0;
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
        if (base_label>=0 && auxsg->seq_vis[v0>>1]!=base_label) {continue;}
        has_overlap = 0;
        if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0, base_label)!=2){
            continue;
        }
        int idx_v = 0;
        uint32_t tmp[2];
        asg_arc_t *tmpav = asg_arc_a(auxsg, v0);
        uint32_t tmpnv = asg_arc_n(auxsg, v0);
        for (uint32_t i=0; i<tmpnv; i++){
            if (tmpav[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[tmpav[i].v>>1]!=base_label) {continue;}
            tmp[idx_v] = tmpav[i].v;
            idx_v++;
        }
        assert(idx_v==2);
        v1 = tmp[0];
        v2 = tmp[1];
        
        if (hamt_asgarc_util_countPre(auxsg, v1, 0, 0, base_label)!=1 || hamt_asgarc_util_countPre(auxsg, v2, 0, 0, base_label)!=1){  // require the bifucation to be simple
            continue;
        }

        // early termination: simple cases
        if (hamt_asgarc_util_checkSimpleBubble(auxsg, v0, NULL, base_label)>0){
            // v0 is the start of a simple bubble
            continue;
        }
        if (hamt_asgarc_util_countNoneTipSuc(auxsg, v1, base_label)<1 || hamt_asgarc_util_countNoneTipSuc(auxsg, v2, base_label)<1){
            // at least one side terminated in a tip
            continue;
        }

        ////// traversal //////
        // (not the most efficient way but this looks cleaner)
        // fprintf(stderr, "traversing utg%.6d dir %d.\n", (int)(v0>>1)+1, (int)(v0&1)); fflush(stderr);
        queue32_reset(&q1); queue32_reset(&q2);
        memset(color1, 0, auxsg->n_seq*2);  memset(color2, 0, auxsg->n_seq*2);
        if (hamt_ugasgarc_util_BFS_mark(ug, v1, color1, &q1, threshold_nread, threshold_bp, base_label)){  // early termination, v1 leads to a not-that-long path
            continue;
        }
        if (hamt_ugasgarc_util_BFS_mark(ug, v2, color2, &q2, threshold_nread, threshold_bp, base_label)){
            continue;
        }
        // (if we reach here, v1 and v2 goes to different & long paths, consider cutting the bifurcation)
        uint32_t san_i;
        for (uint32_t i1=0; i1<auxsg->n_seq*2; i1++){
            // if (color1[i1] && color2[i1]){
            if ( (color1[i1] || color1[i1^1]) && (color2[i1] || color2[i1^1]) ){  // any overlap is not allowed, even it's another direction (which could be the case if there's a bubble being held by a handle, twice)
                if (verbose && (!(color1[i1] && color2[i1]))){
                    fprintf(stderr, "[specialdebug::%s] direction conflict case: utg%.6d utg%.6d utg%.6d\n", __func__,
                                    (int)(v0>>1)+1, (int)(v1>>1)+1, (int)(v2>>1)+1);
                }
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
        fprintf(stderr, "[debug::%s] cut route utg%.6d utg%.6d utg%.6d\n", __func__, (int)(v0>>1)+1, (int)(v1>>1)+1, (int)(v2>>1)+1);
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

int hamt_ugasg_cut_shortTips(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // NOTE: only drop arcs, not deleting sequences
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        if (ug->u.a[vu>>1].len>100000){  // tip not short
            continue;
        }

        // specify the tip direction: no predecessor, has sucessor
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==0 || hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0){
            continue;
        }

        // sancheck + cut: don't drop tip if ALL of its sucessors has no target other than the tip
        //                 i.e. only drop an arc if the target vertex still have other predecessors
        // (note) this intentionally spares all multi-leaf tips; let other rountines try them.
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        int need_to_spare = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            // if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}  // note: do not do this
            wu = av[i].v;
            if (hamt_asgarc_util_countNoneTipSuc(auxsg, wu^1, base_label)==0){
                if (verbose){
                    fprintf(stderr, "[debug::%s] spared an arc for tip utg%.6d (target utg%.6d)\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
                }
            }else{  // cut
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                nb_cut++;
                need_to_spare = 1;
                if (verbose>1){
                    fprintf(stderr, "[debug::%s] cut tip: utg%.6d \n", __func__, (int)(vu>>1)+1);
                }
            }
        }
        // if (!is_hard_drop && !need_to_spare){
        //     hamt_ug_utg_softdel(sg, ug, vu, alt_label);
        // }
    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d unitig tigs.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut;
}

int hamt_ug_cut_shortTips_arbitrary(asg_t *sg, ma_ug_t *ug, int max_length, int base_label){
    // FUNC
    //     cut simple tips shorter than max_length, regardless of whether their targets have no non-tip linkage
    // NOTE
    //     this is meant to be a companion of complex bubble resolve function,
    //     some very short tip may remain after cutting and they may not get treatment
    //     if the complex bubble handle is a tip && is the given tip's target's predecessor.
    //     like
    /*
                     ---(complex bubble edges).........
                    /              .....               \
         start(tip)-------v1 (->complex bubble edges)..end------
                           \
                            ---v0(THE TIP   )
        (since v1^1 apparently has no non-tip targets, v0 will be spared by `hamt_ug_cut_shortTips`,
         which would let this complex bubble remain unresolved in the next round.)
    */
    asg_t *auxsg = ug->g;
    int nb_cut = 0;
    uint32_t nv;
    asg_arc_t *av;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>max_length){continue;}
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==0){continue;}
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            // don't need to check stuff, simply call arc del is fine.
            hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
        }
        nb_cut++;
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d tips (length threshold is %d)\n", __func__, nb_cut, max_length);
    }
    return nb_cut;
}


void hamt_asgarc_markBridges_DFS(asg_t *sg, uint32_t v, uint32_t v_parent, int *DFStime, uint8_t *visited, int *tin, int *low, int base_label){
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
        if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        w = av[i].v;
        if (((v>>1)!=(v_parent>>1)) && (w>>1)==(v_parent>>1)){  // note: v==v_parent is the case for root node
            continue;
        }
        if (visited[w>>1]){
            low[v>>1] = low[v>>1]<tin[w>>1]? low[v>>1] : tin[w>>1];
        }else{
            hamt_asgarc_markBridges_DFS(sg, w, v, DFStime, visited, tin, low, base_label);
            low[v>>1] = low[v>>1]<low[w>>1]? low[v>>1] : low[w>>1];
            if (low[w>>1]>tin[v>>1]){  // arc is bridge
                av[i].is_bridge = 1;
                // mark the complementary arc as bridge too
                sg->arc[asg_arc_get_complementaryID(sg, &av[i])].is_bridge = 1;
            }
        }
    }
}

void hamt_ug_covCutByBridges(asg_t *sg, ma_ug_t *ug, int base_label)
{
    // replacement of the bidirected SCC approach (which did not work as expected)
    // NOTE
    //     the graph is treated as if undirected, aka both directions will be touched when looking for child nodes
    int verbose = 0;
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
            hamt_asgarc_markBridges_DFS(auxsg, v, v, &DFStime, visited, tin, low, base_label);
        }
    }
    if (verbose){
        int n_bridges = 0;
        for (uint32_t i=0; i<auxsg->n_arc; i++){
            if (auxsg->arc[i].is_bridge){
                n_bridges++;
                if (verbose>1){
                    fprintf(stderr, "[debug::%s] identified bridge beween utg%.6d and utg%.6d\n", __func__, (int)(auxsg->arc[i].ul>>33)+1, (int)(auxsg->arc[i].v>>1)+1);
                }
            }
        }
        fprintf(stderr, "[M::%s] %d out of %d arcs are bridge edges.\n", __func__, n_bridges, (int)auxsg->n_arc);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }

    // check arcs
    // (note: since the asg is treated as if undirected, edges in bubbles are not bridges)
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del){continue;}
        if (!auxsg->arc[i].is_bridge){
            continue;
        }
        
        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (base_label>=0 && ( auxsg->seq_vis[v>>1]!=base_label || auxsg->seq_vis[w>>1]!=base_label) ) {continue;}
        if (verbose>1) {fprintf(stderr, "[debug::%s] checking at utg%.6d vs utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);}

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
            if (verbose>1) {fprintf(stderr, "[debug::%s]    failed cov diff: utg%.6d utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);}
            continue;
        }

        // detect circles
        if (hamt_asgarc_detect_circleDFS(auxsg, v, w, 0, base_label)){  // note: v and w is ug with direction
            if (hamt_asgarc_util_isTip(auxsg, v, 0, 0, base_label) && (hamt_asgarc_util_countSuc(auxsg, v, 0, 0, base_label)==1)){
                // don't remove tips from circles
                // (however, allow cutting if the tip is still connected to somewhere else)
                if (verbose>1) {fprintf(stderr, "[debug::%s]    spared bcs handle is a tip\n", __func__);}
                continue;
            }

            // cut
            hamt_ug_arc_del(sg, ug, v, w, 1);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut: utg%.6d - utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                fprintf(stderr, "[debug::%s]     (cov max: %d, cov min:%d)\n", __func__, cov_max, cov_min);
                fprintf(stderr, "[debug::%s]      %d, %d; %d, %d\n", __func__, (int)hamt_asgarc_util_countPre(auxsg, v, 0, 0, base_label), 
                                                                 (int)hamt_asgarc_util_countSuc(auxsg, w, 0, 0, base_label), 
                                                                 (int)hamt_asgarc_util_countPre(auxsg, v, 0, 0, base_label), 
                                                                 (int)hamt_asgarc_util_countSuc(auxsg, w, 0, 0, base_label));
            }
        }else{
            if (verbose) {
                fprintf(stderr, "[debug::%s] spared: utg%.6d - utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                fprintf(stderr, "[debug::%s]    failed circle test\n", __func__);
            }
        }
        
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d links.\n", __func__, nb_cut);
    }
    free(tin);
    free(visited);
    free(low);
    asg_cleanup(sg);
    asg_cleanup(auxsg);
}

int hamt_ug_recover_ovlp_if_existed_core(asg_t *sg, ma_ug_t *ug, uint32_t start_v, uint32_t end_v,
                                    ma_hit_t_alloc *sources, 
                                    const ma_sub_t* coverage_cut, int yes_recover_it){
    // NOTE
    //    start_v and end_v are vertex IDs (with direction)
    //    the corresponding unitig looks like: start_v->.....->end_v (yes, not end_v^1)
    // FUNC
    //    check if overlap ever existed, push or not push it back to the asg
    // RETURN
    //    1 if an arc can be or was recovered 
    //    0 if overlap existed, but the arc is suboptimal (e.g. containment)
    //    -1 if overlap was not available
    int verbose = 0;
    if (verbose){
        fprintf(stderr, "[debug::%s] start %.*s, end %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, start_v>>1), Get_NAME(R_INF, start_v>>1),
                                                                        (int)Get_NAME_LENGTH(R_INF, end_v>>1), Get_NAME(R_INF, end_v>>1));
    }

    int yes_push_arc = 0, yes_push_arc_rev = 0;
    ma_hit_t *h=NULL, *h_rev=NULL;
    if (does_ovlp_ever_exist(sources, end_v, start_v, &h, &h_rev)){
        yes_push_arc = 1;
    }else{
        return -1;
    }
    if (!h || !h_rev) {return 0;}  // (suppress compiler warning)
    asg_arc_t t, *p;
    int ql = coverage_cut[Get_qn(*h)].e - coverage_cut[Get_qn(*h)].s;
    int tl = coverage_cut[Get_tn(*h)].e - coverage_cut[Get_tn(*h)].s;
    int r = ma_hit2arc(h, ql, tl, asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
    if (r>=0){
        // push the 1st direction to the asg
        if (yes_recover_it){
            p = asg_arc_pushp(sg);
            *p = t;
        }
        // push the other direction
        ql = coverage_cut[Get_qn(*h_rev)].e - coverage_cut[Get_qn(*h_rev)].s;
        tl = coverage_cut[Get_tn(*h_rev)].e - coverage_cut[Get_tn(*h_rev)].s;
        r = ma_hit2arc(h_rev, ql, tl, asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t);
        if (r>=0){
            if (yes_recover_it){
                p = asg_arc_pushp(sg);
                *p = t;
            }
            fprintf(stderr, "[debug::%s] success\n", __func__);
            return 1;
        }else{
            fprintf(stderr, "[W::%s] tried to recover arc via a not-so-good overlap (a)\n", __func__);
            return 0;
        }
    }else{
        fprintf(stderr, "[W::%s] tried to recover arc via a not-so-good overlap (b)\n", __func__);
        return 0;
    }    
}


int hamt_ug_recover_ovlp_if_existed(asg_t *sg, ma_ug_t *ug, uint32_t start, uint32_t end,
                                    ma_hit_t_alloc *sources, 
                                    const ma_sub_t* coverage_cut, int search_span){
    // FUNC
    //     check if start->end could've formed an arc. Add the arc accordingly if so
    // PAR
    //     search_span : how many reads to search if the overlap of start/end pair 
    //                     is suboptimal (e.g. containment) & couldn't directly form the arc.
    // NOTE
    //     Topo context of vu will NOT be checked; caller is responsible for sanchecks.
    //     Will search more than start/end vertex - since we can't add arc if the overlap was containment etc.
    //       if all recoverble overlaps were not ideal, do nothing.
    //     After recovering the arc, also trim off reads involved in containment. (reads will be DELETED)
    // RETURN
    //     1 upon recovereing an arc (+cleanup)
    //     0 if no termination criteria met, but ran out of search_span
    //     -1 if we shouldn't recover stuff
    int verbose = 1;

    uint32_t v, w, v_start, v_end;
    asg_t *auxsg = ug->g;
    int status;

    v = start;
    w = end;
    status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, 1);
    if (status<0){
        return -1;
    }else if (status>0){
        fprintf(stderr, "ideal recovery\n");
        return 1;
    }
    if (search_span<0){  // for more generic usage
        return 0;
    }

    // there was overlap, but wasn't able to form the arc
    // search other pairs: step at the unitig end
    if (hamt_asgarc_util_countSuc(sg, v^1, 0, 0, -1)!=1){return 0;}  // assuming v is the end^1 vertex of a unitig (w is a unitig start)
    if (hamt_asgarc_util_get_the_one_target(sg, v^1, &v, 0, 0, -1)<0) {
        fprintf(stderr, "[W::%s] hamt_asgarc_util_get_the_one_target failed\n", __func__);
        return -1;
    }
    v ^= 1;
    for (int i=0; i<search_span-1; i++){
        status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, 0);
        if (status>0){
            // step the other side
            uint32_t w_end = w;
            int idx = 0;
            for (int j=0; j<search_span-1; j++){
                // step
                if (hamt_asgarc_util_countSuc(sg, w, 0, 0, -1)!=1){break;}
                if (hamt_asgarc_util_get_the_one_target(sg, w, &w, 0, 0, -1)<0){break;}
                // check
                status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, 0);
                if (status>0){
                    w_end = w;
                    idx = j;
                }else if (status<0){break;}
            }
            if (idx==(search_span-1)){fprintf(stderr, "[W::%s] search span not big enough?\n", __func__);}
            status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w_end, sources, coverage_cut, 1);
            assert(status==1);
            fprintf(stderr, "stepped recovery, read pair: %.*s - %.*s\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1),
                                                                           (int)Get_NAME_LENGTH(R_INF, w_end>>1), Get_NAME(R_INF, w_end>>1));
            return 1;
        }else if (status<0){  // encounter any non-overlapping pairs and we're done
            break;
        }
        if (hamt_asgarc_util_countSuc(sg, v^1, 0, 0, -1)!=1){return 0;}
        if (hamt_asgarc_util_get_the_one_target(sg, v^1, &v, 0, 0, -1)<0) {
            fprintf(stderr, "[W::%s] hamt_asgarc_util_get_the_one_target failed\n", __func__);
            return -1;
        }
        v ^= 1;
    }

    v = start;
    w = end;
    // search other pairs: step at the unitig start
    if (hamt_asgarc_util_countSuc(sg, w, 0, 0, -1)!=1){return 0;}  // assuming v is the end^1 vertex of a unitig (w is a unitig start)
    if (hamt_asgarc_util_get_the_one_target(sg, w, &w, 0, 0, -1)<0) {
        fprintf(stderr, "[W::%s] hamt_asgarc_util_get_the_one_target failed\n", __func__);
        return -1;
    }
    for (int i=0; i<search_span-1; i++){
        status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, 1);
        if (status>0){
            // step the other side
            uint32_t v_end = v;
            int idx = 0;
            for (int j=0; j<search_span-1; j++){
                // step
                if (hamt_asgarc_util_countSuc(sg, v^1, 0, 0, -1)!=1){break;}
                if (hamt_asgarc_util_get_the_one_target(sg, v^1, &v, 0, 0, -1)<0){break;}
                v^=1;
                // check
                status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, 0);
                if (status>0){
                    v_end = v;
                    idx = j;
                }else if (status<0){break;}
            }
            if (idx==(search_span-1)){fprintf(stderr, "[W::%s] search span not big enough?\n", __func__);}
            status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v_end, w, sources, coverage_cut, 1);
            assert(status==1);
            fprintf(stderr, "stepped recovery, read pair: %.*s - %.*s\n", (int)Get_NAME_LENGTH(R_INF, v_end>>1), Get_NAME(R_INF, v_end>>1),
                                                                           (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1));
            return 1;
        }else if (status<0){  // encounter any non-overlapping pairs and we're done
            return -1;
        }
        if (hamt_asgarc_util_countSuc(sg, w, 0, 0, -1)!=1){return 0;}
        if (hamt_asgarc_util_get_the_one_target(sg, w, &w, 0, 0, -1)<0) {
            fprintf(stderr, "[W::%s] hamt_asgarc_util_get_the_one_target failed\n", __func__);
            return -1;
        }
    }
    return 0;
}


//////////////////////////////////////////////////////////////
//               resolve tangles etc                        //
//////////////////////////////////////////////////////////////
int hamt_ug_check_complexBubble(asg_t *sg, ma_ug_t *ug, int max_size, uint32_t v0, uint32_t *end0, int base_label){
    // FUNC
    //    use a BFS variant to test if v0 is the start of a complex bubble (of size less than `max_size`)
    //    (complex bubble means: no inversion aka all vertices shall be discovered in only 1 direction,
    //                          no loop/circle)
    //    if so, arbitrarily chose a path from it and drop all the other arcs.
    // RETURN
    //    1 if yes (side effect: store the end vertex in *end0)
    //    0 if no
    // NOTE
    //    not sure if this is always correct. Haven't seen it break in practice, but have no proof.

    int verbose = 2;

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, v_tmp, nv, nw, tmp_nv;
    asg_arc_t *av, *aw, *tmp_av;
    int nb_nodes_visited = 0;
    vu = v0;

    if (verbose) {fprintf(stderr, "[debug::%s] check vu utg%.6d\n", __func__, (int)(v0>>1)+1);}

    // basic topo checks, requiring vu to be a articulation-vertex-like vertex
    // (note that to this point, we should have no trivial tips.)
    // (forward direction)
    if (hamt_asgarc_util_countNoneTipSuc(auxsg, vu, base_label)<2){return 0;}
    nv = asg_arc_n(auxsg, vu);
    av = asg_arc_a(auxsg, vu);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)!=1){
            // allow backward tip
            int san = 0;
            tmp_nv = asg_arc_n(auxsg, av[i].v^1);
            tmp_av = asg_arc_a(auxsg, av[i].v^1);
            for (int i_tmp=0; i_tmp<tmp_nv; i_tmp++){
                if (tmp_av[i_tmp].del){continue;}
                if (base_label>=0 && auxsg->seq_vis[tmp_av[i_tmp].v>>1]!=base_label) {continue;}
                if (tmp_av[i_tmp].v^1 == vu){continue;}  // note: must check direction
                if (hamt_asgarc_util_countSuc(auxsg, tmp_av[i_tmp].v, 0, 0, base_label)>0){  // the backward bifurcation leads to non-tip tig
                    san=1;
                    break;
                }
            }
            if (san){
                if (verbose) {fprintf(stderr, "[debug::%s]     exit bc target backward bifur\n", __func__);}
                return 0;
            }
        }
    }
    // // (backward direction)
    // if (hamt_asgarc_util_countPre(auxsg, vu^1, 0, 0)>0){
    //     nv = asg_arc_n(auxsg, vu^1);
    //     av = asg_arc_a(auxsg, vu^1);
    //     for (uint32_t i=0; i<nv; i++){
    //         if (av[i].del){continue;}
    //         v_tmp = av[i].v^1;
    //         if (hamt_asgarc_util_countSuc(auxsg, v_tmp, 0, 0)!=1){
    //             return 0;
    //         }
    //     }
    // }

    // BFS variant init
    queue32_t q[2];
    queue32_init(&q[0]);
    queue32_init(&q[1]);
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    queue32_enqueue(&q[0], vu);
    color[v0] = 1;
    nb_nodes_visited++;
    uint8_t iq = 0;  // which queue is the "working" one

    // (using the stack as a simple vector to store current round's color info)
    stacku32_t stack;
    stacku32_init(&stack);

    int is_abnormal = 0, is_tangle = 0, nb_rounds = 0;
    int has_anything_happend = 1;
    int v0_special_flag = 1;
    

    // test if we have a very simple tangle here
    // a BFS variant:
    //     We have 2 queues. For each turn, check every element in queue a,
    //     for each white child node, if all of the child node's parents is seen *before the current round*,
    //     (tmporarily store the current round's discoveries in a seperate vector color_tmp)
    //     dequeue it like normal BFS; otherwise, enqueue it to queue b. For the next turn,
    //     the current queue b will be queue a, and current a will be b.
    //     (apply color_tmp to color)
    // special note
    //     if the walk loops, only qualify the bubble if handle is long enough and the bubble is short 
    //     (otherwise it's the broken spoon subgraph in mock)
    //     Limited by the current test data, we just check the handle length and let 1Mbp or something be the threshold
    //     TODO - improve it.
    int longest_length_in_bubble = 0;
    while (1){
        // termination criteria
        if (!has_anything_happend){
            if (verbose) {fprintf(stderr, "[debug::%s]     exiting bc nothing happend last round\n", __func__);}
            break;
        }
        if (is_abnormal){
            if (verbose) {fprintf(stderr, "[debug::%s]     exiting bc abnormal\n", __func__);}
            break;
        }
        if (queue32_get_size(&q[iq])==1){  // found the end vertex
            if (nb_nodes_visited>1){
                /////// update end0 ///////
                queue32_dequeue(&q[iq], end0);
                is_tangle = 1;
                if (verbose) {fprintf(stderr, "[debug::%s]     found end vertex utg%.6d\n", __func__, (int)((*end0)>>1)+1);}
                break;
            }
        }
        if (queue32_get_size(&q[iq])==0){  // queue emptied without finding the end vertex
            if (verbose) {fprintf(stderr, "[debug::%s]     exiting bc queue empty\n", __func__);}
            break;
        }
        if (nb_nodes_visited>max_size){  // time out
            if (verbose) {fprintf(stderr, "[debug::%s]     exiting bc exceeding max nodes\n", __func__);}
            break;
        }
        
        // current round
        has_anything_happend = 0;
        nb_rounds++;
        if (nb_rounds>1000){
            fprintf(stderr, "debug break: more than 1000 round\n");
            exit(1);
        }
        stacku32_reset(&stack);  // holds color info of the current round
        int vu_not_passing;
        if (verbose) {fprintf(stderr, "[debug::%s]     > a new round\n", __func__);}
        while (queue32_dequeue(&q[iq], &vu)){
            if (verbose) {fprintf(stderr, "[debug::%s]     at utg%.6d\n", __func__, (int)(vu>>1)+1);}
            if (color[vu]==2){continue;}
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);

            vu_not_passing = 0;
            // 1st pass: check if all predecessors of wu have been seen. If any wu not passing, put vu back.
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del) {continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                wu = av[i].v;
                if (color[wu]!=0) {continue;}
                nw = asg_arc_n(auxsg, wu^1);
                aw = asg_arc_a(auxsg, wu^1);
                uint32_t iw;
                for (iw=0; iw<nw; iw++){
                    if (aw[iw].del){continue;}
                    if (base_label>=0 && auxsg->seq_vis[aw[iw].v>>1]!=base_label) {continue;}
                    if (hamt_asgarc_util_countSuc(auxsg, aw[iw].v, 0, 0,base_label)==0 &&
                        hamt_asgarc_util_countPre(auxsg, aw[iw].v, 0, 0, base_label)==1) {  // ignore simple tip predecessors
                        color[aw[iw].v] = 2;
                    }else{
                        if (!(color[aw[iw].v] || color[aw[iw].v^1])){  // not all predecessors have been seen; take back vu (if not already have) and don't push wu
                            vu_not_passing = 1;
                            if (verbose) {fprintf(stderr, "[debug::%s]     -> wu utg%.6d did not pass pred check\n", __func__, (int)(wu>>1)+1);}
                            break;
                        }
                    }
                }
                if (vu_not_passing){break;}
            }
            if (vu_not_passing){
                queue32_enqueue(&q[iq^1], vu);  // take back vu
                if (verbose) {fprintf(stderr, "[debug::%s]     take back vu\n", __func__);}
            }else{
                stacku32_push(&stack, vu);  // mark color (merge to `color` later)
                has_anything_happend = 1;  // since vu is finished
            }
            // 2nd pass: push passing white child node
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del) {continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                wu = av[i].v;
                if (color[wu]!=0) {continue;}
                // set to finished if the child is a simple tip (in the current walking direction)
                if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)==0 &&
                    hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)==1){  // note: must check pre, since the actual end vertex might be a tip (with more than 1 pre)
                    color[wu] = 2;
                    continue;
                }
                nw = asg_arc_n(auxsg, wu^1);
                aw = asg_arc_a(auxsg, wu^1);
                uint32_t iw;
                for (iw=0; iw<nw; iw++){
                    if (aw[iw].del){continue;}
                    if (base_label>=0 && auxsg->seq_vis[aw[iw].v>>1]!=base_label) {continue;}
                    // on pardon predecessor-unseen-but-its-a-tip condition: it should already have been colored black
                    if (!(color[aw[iw].v] || color[aw[iw].v^1])){ // predecessor check
                        break;
                    }
                }
                if (iw==nw){  // passed predecessor checking
                    // check if wu has already been pushed in this round
                    if (stack32_is_in_stack(&stack, wu)){
                        continue;
                    }
                    // push
                    queue32_enqueue(&q[iq^1], wu);
                    if (((wu>>1)!=(vu>>1)) && ((wu>>1)!=(v0>>1)) && ug->u.a[wu>>1].len>longest_length_in_bubble) {   // log length info
                        longest_length_in_bubble = ug->u.a[wu>>1].len ; 
                    }
                    stacku32_push(&stack, wu);
                    nb_nodes_visited++;
                    has_anything_happend = 1;  // since we pushed a new node
                    if (verbose) {fprintf(stderr, "[debug::%s]     push utg%.6d\n", __func__, (int)(wu>>1)+1);}
                }
            }
        }
        // update color
        while (stacku32_pop(&stack, &vu)){
            if (vu==v0 && v0_special_flag){
                // A special case: after v0 is marked finished for the 1st time (the initial push), reset its color to white. 
                // This is supposed to enable detecting complex bubble in a circle (i.e. start and end is both v0)
                v0_special_flag = 0;
                color[v0] = 0;
                continue;
            }
            color[vu]++;
            if (color[vu^1]){  // sancheck: only one direction shall be discovered
                is_abnormal = 1;
                break;
            }
        }
        stacku32_reset(&stack);

        // switch working queue
        queue32_reset(&q[iq]);
        iq ^= 1;
    }

    // check length 
    if (is_tangle && (v0==(*end0))){
        if (ug->u.a[v0>>1].len<(longest_length_in_bubble*2)){
            is_tangle = 0;
            if (verbose){
                fprintf(stderr, "[debug::%s] detected complex bubble for handle utg%.6d but discarded bc length (v0 %d, longest %d)\n", __func__, (int)(v0>>1)+1, (int)ug->u.a[v0>>1].len, longest_length_in_bubble);
            }
        }
    }

    ////////// debug sanchecks /////////////
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){  // color shouldn't be larger than 2
        if (color[i]>2){
            fprintf(stderr, "[debugERROR::%s] color can't be larger than 2 (utg%.6d)\n", __func__, (int)(i>>1)+1);fflush(stderr);
            exit(1);
        }
    }
    ////////////////////////////////////////

    queue32_destroy(&q[0]);
    queue32_destroy(&q[1]);
    stacku32_destroy(&stack);
    free(color);

    return is_tangle;
}

int hamt_ug_pop_complexBubble(asg_t *sg, ma_ug_t *ug, uint32_t start0, uint32_t end0, 
                              int base_label, int alt_label, int is_hard_drop){
    // FUNC
    //     start0 and end0 is the start/end handle of a complex bubble (per `hamt_ug_check_complexBubble`)
    //     retain an arbitrary path (greedily long but not the longest) and discard other arcs
    // TODO
    //     try to be not too arbitrary
    // TODO 2
    //    maybe instead of dropping arcs, i should remove the seqs?
    //    (does hamt need alternative ctg? yes)
    // RETURN
    //     number of dropped arcs
    int verbose = 2;

    asg_t *auxsg = ug->g;
    int nb_cut = 0;
    uint32_t vu, wu, nv, idx;
    asg_arc_t *av;
    vu = start0;
    int has_anything_happened = 1;

    // DFS
    stacku32_t stack;
    stacku32_init(&stack);
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    color[start0] = 1;
    int is_exhausted = 0, yes_break = 0;

    // for backtrace
    stacku32_t stack_bt;
    stacku32_init(&stack_bt);

    // for length selection
    uint64_t buf[50];
    uint32_t node_length;

    // init
    stacku32_push(&stack, start0);
    stacku32_push(&stack_bt, start0);

    while (stacku32_pop(&stack, &vu)){
        if (verbose>1){fprintf(stderr, "        pop, at utg%.6d\n", (int)(vu>>1)+1);}
        if (color[vu]==2) {continue;}
        nv = asg_arc_n(auxsg, vu); 
        av = asg_arc_a(auxsg, vu);

        // is all child exhausted?
        is_exhausted = 1;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (av[i].v==end0){  // 
                yes_break = 1;
                stacku32_push(&stack_bt, end0);
                vu = end0;
                if (verbose){fprintf(stderr, "[debug::%s] reached the end properly(1)\n", __func__);}
                break;
            }
            if (color[av[i].v]==0){
                // if the unseen child is a tip (in the current walking direction), set it to finished
                if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label)==0 &&
                    hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)==1){  // need to count predecessor, in case the end happens to be a tip
                    color[av[i].v] = 2;
                }else{
                    is_exhausted = 0;
                    break;
                }
            }
        }
        if (yes_break) {break;}

        if (is_exhausted){
            if (verbose>1){fprintf(stderr, "           exhausted utg%.6d\n", (int)(vu>>1)+1);}
            if (!stacku32_pop(&stack_bt, &wu)){  // remove the current node from the recorded path
                fprintf(stderr, "[W::%s] recorded path went empty, abort\n", __func__);
                stacku32_destroy(&stack);
                stacku32_destroy(&stack_bt);
                free(color);
                return nb_cut;
            }
            color[vu] = 2;
            continue;
        }else{
            // put back the current node
            stacku32_push(&stack, vu);
            if (verbose>1){fprintf(stderr, "           put back utg%.6d\n", (int)(vu>>1)+1);}
            // pick the longest child
            if (nv>50){  // boundary check
                fprintf(stderr, "[W::%s] aborted cutting because number child nodes is larger than static buffer.\n", __func__);
                break;
            }
            int idx = 0;
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del) {continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                if (color[av[i].v]!=0) {continue;}
                if (av[i].v==end0){  // found the end node, terminate
                    stacku32_push(&stack_bt, end0);
                    vu = end0;
                    if (verbose){fprintf(stderr, "[debug::%s] reached the end properly(2)\n", __func__);}
                    yes_break = 1;
                    break;
                }
                node_length = ug->u.a[av[i].v>>1].len;
                buf[idx++] = ((uint64_t)node_length)<<32 | ((uint64_t)av[i].v);
            }
            if (yes_break) {break;}  // found the end node
            assert(idx);
            // push
            radix_sort_ovhamt64(buf, buf+idx);
            stacku32_push(&stack_bt, (uint32_t)buf[idx-1]);
            for (int i=0; i<idx; i++){  // push shorter child node first
                stacku32_push(&stack, (uint32_t)buf[i]);
                color[(uint32_t)buf[i]] = 1;
                if (verbose>1){fprintf(stderr, "           push utg%.6d\n", (int)(((uint32_t)buf[i])>>1)+1);}
            }
        }
    }
    if (vu!=end0){  // sancheck, shouldn't happen
        fprintf(stderr, "[E::%s] weren't able to reach the end of the complex bubble. Traceback:\n", __func__);
        fprintf(stderr, "[E::%s] start utg%.6d , end utg%.6d ; recorded path (backward):\n", __func__, (int)(start0>>1)+1, (int)(end0>>1)+1);
        while (stacku32_pop(&stack_bt, &vu)){
            fprintf(stderr, "[E::%s]     utg%.6d\n", __func__, (int)(vu>>1)+1);
        }
        fprintf(stderr, "[E::%s] (end of traceback)\n", __func__);
        stacku32_destroy(&stack);
        stacku32_destroy(&stack_bt);
        free(color);
        return 0;
    }

    // drop all other arcs
    if (is_hard_drop){
        nb_cut = hamt_ug_arc_del_between(sg, ug, start0, end0, 1, base_label);
    }else{
        nb_cut = hamt_ug_arc_flip_between(sg, ug, start0, end0, base_label, alt_label);
    }
    // recover the stored path
    assert(stacku32_pop(&stack_bt, &vu));
    if (verbose){
        fprintf(stderr, "[debug::%s] retaining the following path between utg%.6d and utg%.6d\n", __func__, (int)(start0>>1)+1, (int)(end0>>1)+1);
    }
    while (stacku32_pop(&stack_bt, &wu)){
        if (verbose){
            fprintf(stderr, "[debug::%s]     utg%.6d -> utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
        }
        if (is_hard_drop){
            hamt_ug_arc_del(sg, ug, wu, vu, 0);
        }else{
            hamt_ug_arc_del(sg, ug, wu, vu, 0);
            hamt_ug_utg_softdel(sg, ug, wu, base_label);
            hamt_ug_utg_softdel(sg, ug, vu, base_label);
        }
        vu = wu;
        nb_cut--;
    }
    // hamt_ug_cleanup_arc_by_labels(sg, ug);

    stacku32_destroy(&stack);
    stacku32_destroy(&stack_bt);
    free(color);
    return nb_cut;
}


//////////////////////////////////////////////////////////
//                   interface                          //
//////////////////////////////////////////////////////////


void hamt_circle_cleaning(asg_t *sg, ma_ug_t *ug, int base_label)
{
    // hamt_asgarc_ugCovCutSCC(sg, ug, coverage_cut, sources, ruIndex);  
    hamt_ug_covCutByBridges(sg, ug, base_label);
    hamt_asgarc_ugCovCutDFSCircle(sg, ug, base_label);
}

void hamt_clean_shared_seq(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    hamt_asgarc_ugTreatMultiLeaf(sg, ug, 50000, base_label, alt_label, is_hard_drop);
    hamt_asgarc_ugTreatBifurcation(sg, ug, -1, 10000000, base_label);
}

int hamt_ug_pop_bubble(asg_t *sg,ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop)
{
    // remove simple bubbles at unitig graph level
    // NOTE: only drop arcs, not deleting sequences
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu[2], uu, nv;
    asg_arc_t *av;
    int cov[2], idx, nb_cut=0, nb_tip_cut;

    // pre clean
    nb_tip_cut = hamt_ugasg_cut_shortTips(sg, ug, base_label, alt_label, is_hard_drop);  // note: does not affect multi-leaf tips

    // simple bubble popping
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=2){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleBubble(auxsg, vu, &uu, base_label)==0){  // check if the current utg is the start of a simple bubble
            continue;
        }
        // collect the bubble edges
        av = asg_arc_a(auxsg, vu);
        for (int i=0, i_=0; i<asg_arc_n(auxsg, vu); i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            wu[i_++] = av[i].v;
            assert(i_<=2);
        }
        
        // pop
        if (ug->u.a[wu[0]>>1].len>100000 || ug->u.a[wu[1]>>1].len>100000){
            // if either side of the bubble is long,
            // pop the shorter side
            idx = ug->u.a[wu[0]>>1].len > ug->u.a[wu[1]>>1].len ? 1 : 0;
        }else{
            // otherwise, pop the side with lesser coverage
            idx = ug->utg_coverage[wu[0]>>1] > ug->utg_coverage[wu[1]>>1]? 1:0;
        }
        
        // drop arc
        if (verbose>1){
            fprintf(stderr, "[debug::%s]> utg%.6d\n", __func__, (int)(vu>>1)+1);
            fprintf(stderr, "             vu is utg%.6d, wu is %.6d, uu is %.6d\n", (int)(vu>>1)+1, (int)(wu[idx]>>1)+1, (int)(uu>>1)+1);
        }
        if (is_hard_drop){
            hamt_ug_arc_del(sg, ug, vu, wu[idx], 1);
            hamt_ug_arc_del(sg, ug, wu[idx], uu, 1);
        }else{
            hamt_ug_arc_del(sg, ug, vu, wu[idx], 1);
            hamt_ug_arc_del(sg, ug, wu[idx], uu, 1);
            hamt_ug_utg_softdel(sg, ug, wu[idx], alt_label);
        }
        nb_cut++;
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        asg_symm(sg);
        asg_symm(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d simple bubbles on unitig graph.\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut+nb_tip_cut;
}

int hamt_ug_pop_miscbubble(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //    pop small multi-edge bubble, or bubbles cluttered together,  tolerating tips
    // NOTE
    //     affected topo should include:
    //       - multi-edge "bubbles"
    //       - not-strictly-a-bubble bubbles (i.e. start/end node has other linkage)
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu[20], uu[20];
    uint32_t nv;
    int idx;
    int color[20], flag, nb_cut = 0;
    asg_arc_t *av;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        // // spare larger structures
        // if (ug->u.a[vu>>1].len>50000){
        //     continue;
        // }

        // specify direction: vu shall be the start of the misc bubble
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<2){
            continue;
        }
        idx = 0;
        memset(color, 0, sizeof(int)*20);
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label)!=1 || hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)!=1){
                continue;
            }
            // if (ug->u.a[av[i].v>>1].len>50000){
            //     continue;
            // }
            wu[idx] = av[i].v;
            hamt_asgarc_util_get_the_one_target(auxsg, av[i].v, &uu[idx], 0, 0, base_label);
            idx++;
            if (idx==20){  // buffer sancheck
                fprintf(stderr, "[W::%s] more than 20 targets?? \n", __func__);
                break;
            }
        }
        if (idx<2){
            continue;
        }

        // check if any edge shares uu with another edge; preserve the 1st one and discard all others
        // TODO: coverage-aware
        if (verbose){
            fprintf(stderr, "[debug::%s] check at utg%.6d\n", __func__, (int)(vu>>1)+1);
        }
        int n_candidates = 1;
        uint32_t buf_candidates[20];
        int max_l = 0, idx_candidate;
        int i1;
        for (i1=0; i1<idx; i1++){
            if (n_candidates>1){
                max_l = 0;
                idx_candidate = 0;
                for (int ic=0; ic<n_candidates; ic++){
                    if (ug->u.a[buf_candidates[ic]>>1].len>max_l){
                        max_l = ug->u.a[buf_candidates[ic]>>1].len;
                        idx_candidate = ic;
                    }
                }
                for (int ic=0; ic<n_candidates; ic++){
                    if (ic==idx_candidate){continue;}
                    hamt_ug_arc_del(sg, ug, vu, buf_candidates[ic], 1);
                    hamt_ug_arc_del(sg, ug, buf_candidates[ic], uu[i1-1], 1);
                    if (verbose){
                    fprintf(stderr, "[debug::%s] cut, u utg%.6d  w utg%.6d u utg%.6d\n", __func__,
                                        (int)(vu>>1)+1, (int)(buf_candidates[ic]>>1)+1, (int)(uu[i1-1]>>1)+1);
                    }
                    nb_cut++;
                }
            }
            n_candidates = 1;
            buf_candidates[0] = wu[i1];
            if (color[i1]){continue;}  // this wu has been examined/cut
            // if (hamt_asgarc_util_isTip(auxsg, uu[i1], 0, 0)){continue;}  // ignore misc bubble that ends at a tig

            for (int i2=i1+1; i2<idx; i2++){
                if (color[i2]){continue;}
                if (uu[i2]!=uu[i1]){continue;}
                color[i2] = 1;
                buf_candidates[n_candidates] = wu[i2];
                n_candidates++;
            }
        }
    }
    // cleanup
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d spots\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut;
}

int hamt_ug_pop_miscbubble_aggressive(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //    pop small multi-edge bubble, or bubbles cluttered together,  
    //    not only tolerating tips, but also branchings
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, wu2;
    uint32_t nv, nw, nw2;
    asg_arc_t *av, *aw, *aw2;
    int flag = 0, nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        // spare larger structures
        if (ug->u.a[vu>>1].len>50000){
            continue;
        }
        // specify direction: vu shall be the start of the misc bubble
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<2){
            continue;
        }
        
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<(nv-1); i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label)==0){continue;}  // ignore tip
            if (ug->u.a[av[i].v>>1].len>50000){continue;}  // ignore larger structures
            wu = av[i].v;

            nw = asg_arc_n(auxsg, wu);
            aw = asg_arc_a(auxsg, wu);
            for (uint32_t iw=0; iw<nw; iw++){  // check wu's targets
                for (uint32_t i2=i+1; i2<nv; i2++){  // check other targets (aka wu2) of vu, and see if any of them will end up sharing a target with wu
                    if (av[i2].del){continue;}
                    if (base_label>=0 && auxsg->seq_vis[av[i2].v>>1]!=base_label) {continue;}
                    if (hamt_asgarc_util_countSuc(auxsg, av[i2].v, 0, 0, base_label)==0){continue;}  // ignore tip
                    if (ug->u.a[av[i2].v>>1].len>50000){continue;}  // ignore larger structures

                    wu2 = av[i2].v;
                    nw2 = asg_arc_n(auxsg, wu2);
                    aw2 = asg_arc_a(auxsg, wu2);
                    for (uint32_t iw2=0; iw2<nw2; iw2++){  // check the other target's target
                        if (aw2[iw2].del){continue;}
                        if (base_label>=0 && auxsg->seq_vis[aw2[iw2].v>>1]!=base_label) {continue;}
                        if (aw2[iw2].v!=aw[iw].v){continue;}
                        // cut
                        // note: don't cut between vu and wu2
                        hamt_ug_arc_del(sg, ug, wu2, aw[iw].v, 1);
                        if (verbose){
                            fprintf(stderr, "[debug::%s] cut, u utg%.6d w utg%.6d u utg%.6d\n", __func__,
                                    (int)(vu>>1)+1, (int)(wu2>>1)+1, (int)(aw[iw].v>>1)+1);
                        }
                        nb_cut++;
                    }

                }

            }

        }
        
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d spots\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut;
}


int hamt_ug_pop_simpleInvertBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // TODO: this function isn't working properly.
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu;
    int nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>50000){
            continue;
        }
        if (hamt_asgarc_util_checkSimpleInvertBubble(sg, vu, &wu, &uu, base_label)){  // note: uu's direction is adjusted for popping; no need to ^1
            if (ug->u.a[wu>>1].len>50000){
                continue;
            }
            if (ug->u.a[uu>>1].len>50000){
                continue;
            }
            if (is_hard_drop){
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                hamt_ug_arc_del(sg, ug, wu, uu, 1);
            }else{
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                hamt_ug_arc_del(sg, ug, wu, uu, 1);
                hamt_ug_utg_softdel(sg, ug, wu, alt_label);
            }
            nb_cut++;

            if (verbose){
                fprintf(stderr, "[debug::%s] v utg%.6d w utg%.6d u utg%.6d \n", __func__, 
                                (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
            }
        }
    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d locations\n", __func__, nb_cut);
        fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    return nb_cut;
}



int hamt_ug_pop_tinyUnevenCircle(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // experimental
    /*
             -----v1-----
           /     /  \    \
    -----v0     |   |     v3-----
                \  /
                 v2
        if v0, v1 and v3 have similar ok-ish coverages, while v2's is very low,
        let's drop v2, it's not very likely to be a repeat element.
    */
   // NOTE
   //    coverage criteria used are defined in this function.
   int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t v[4];
    uint32_t w1, w2, nv;
    asg_arc_t *av;
    int cov[4], diff0, diff1, diff2;
    int nb_cut = 0;
    int san_flag;

    for (v[2]=0; v[2]<auxsg->n_seq*2; v[2]++){
        if (base_label>=0 && auxsg->seq_vis[v[2]>>1]!=base_label){continue;}

        // shall be short
        if (ug->u.a[v[2]>>1].len>50000){continue;}
        if (verbose>1){fprintf(stderr, "[debug::%s] at utg%.6d\n", __func__, (int)(v[2]>>1)+1);}

        // topo of v2
        if (hamt_asgarc_util_countPre(auxsg, v[2], 0, 0, base_label)!=1 || hamt_asgarc_util_countSuc(auxsg, v[2], 0, 0, base_label)!=1){
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed 1\n", __func__);}
            continue;
        }   

        // try to get v1
        if (hamt_asgarc_util_get_the_one_target(auxsg, v[2], &w1, 0, 0, base_label)<0){fprintf(stderr, "ERROR %s\n", __func__); fflush(stderr); exit(1);}
        if (hamt_asgarc_util_get_the_one_target(auxsg, v[2]^1, &w2, 0, 0, base_label)<0){fprintf(stderr, "ERROR %s\n", __func__); fflush(stderr); exit(1);}
        // topo check if v1 is a such v1, and whether v0 and v3 both check
        if (w1!=(w2^1)){  // check v1
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed 2\n", __func__);}
            continue;
        }
        if ((w1>>1)==(v[2]>>1)){  // self loop
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed 3\n", __func__);}
            continue;
        }
        if (hamt_asgarc_util_countSuc(auxsg, w1, 0, 0, base_label)!=2 || hamt_asgarc_util_countPre(auxsg, w1, 0, 0, base_label)!=2){
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed 4\n", __func__);}
            continue;
        }
        v[1] = w1;

        // get v0 and v3
        nv = asg_arc_n(auxsg, v[1]);
        av = asg_arc_a(auxsg, v[1]);
        san_flag = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            if ((av[i].v>>1)==(v[2]>>1)) {continue;}
            v[0] = av[i].v;
            san_flag = 1;
            break;
        }
        if (!san_flag) {
            fprintf(stderr, "[E::%s] failed to get v0 when supposed to\n", __func__);
            continue;
        }
        nv = asg_arc_n(auxsg, v[1]^1);
        av = asg_arc_a(auxsg, v[1]^1);
        san_flag = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            if ((av[i].v>>1)==(v[2]>>1)) {continue;}
            v[3] = av[i].v;
            san_flag = 1;
            break;
        }
        if (!san_flag) {
            fprintf(stderr, "[E::%s] failed to get v3 when supposed to\n", __func__);
            continue;
        }
        
        // sancheck: v0 and v3 shall not be v1 or v2. 
        if ((v[0]>>1)==(v[1]>>1) || (v[0]>>1)==(v[2]>>1) || (v[3]>>1)==(v[1]>>1) || (v[3]>>1)==(v[2]>>1)){  // compiler note: is safe
            if (verbose>1) {fprintf(stderr, "[debug::%s]     failed 5\n", __func__);}
            continue;
        }

        // get coverages
        for (int i=0; i<4; i++){
            cov[i] = ug->utg_coverage[v[i]>>1];
        }
        // check coverages
        if (verbose>1){
            fprintf(stderr, "[debug::%s]     check coverage\n", __func__);
            fprintf(stderr, "[debug::%s]       coverages were: %d, %d, %d; %d\n", __func__, cov[0], cov[1], cov[3], cov[2]);
            fprintf(stderr, "[debug::%s]       vertices were: %d, %d, %d; %d\n", __func__, v[0], v[1], v[3], v[2]);
        }
        if ( ((cov[1]-cov[2])>=5) && (((float)cov[1]/cov[2]>=2) || cov[2]==0)){  // note: ha's get utg coverage could give 0, so waterproofing division by zero here
            if (verbose>1){fprintf(stderr, "[debug::%s]       cov checkpoint 1\n", __func__);}
            if ( (cov[0]>cov[2]) && (cov[3]>cov[2]) ) {
                if (verbose>1){fprintf(stderr, "[debug::%s]       cov passed\n", __func__);}
                // cut
                if (is_hard_drop){
                    hamt_ug_arc_del(sg, ug, v[1], v[2], 1);
                    hamt_ug_arc_del(sg, ug, v[1]^1, v[2]^1, 1);
                }else{
                    hamt_ug_arc_del(sg, ug, v[1], v[2], 1);
                    hamt_ug_arc_del(sg, ug, v[1]^1, v[2]^1, 1);
                    hamt_ug_utg_softdel(sg, ug, v[2], alt_label);
                }
                if (verbose){
                    fprintf(stderr, "[debug::%s] at utg%.6d , removed its utg%.6d\n", __func__, (int)(v[1]>>1)+1, (int)(v[2]>>1)+1);
                    fprintf(stderr, "[debug::%s]     coverages were: %d, %d, %d; %d\n", __func__, cov[0], cov[1], cov[3], cov[2]);
                }
                nb_cut++;
            }
        }
   }
   if (VERBOSE){
       fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
   }
   return nb_cut;

}



int hamt_ug_pop_terminalSmallTip(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // FUNC
    //    if all targets of a unitig are tips (*tips that's only connect to this one unitig!)
    //    keep only the longest one
    // NOTE
    //    doesn't care about haplotype, just arbitrarily take the longest one.
    //    Name says "small" but isn't doing length check; in meta they just appeared to be all rather short.
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv;
    asg_arc_t *av;
    int i, idx, l;
    int nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu); 
        l = 0;
        for (i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            wu = av[i].v;
            if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)!=0 || hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)!=1){
                break;
            }
            if (ug->u.a[wu>>1].len>=l){
                l = ug->u.a[wu>>1].len;
                idx = i;
            }
        }
        if (i!=nv){continue;}  // not all targets are tips
        for (i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (i==idx) {continue;}  // compile note: is safe  // keep the longest tip
            wu = av[i].v;
            if (is_hard_drop){
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
            }else{
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                hamt_ug_utg_softdel(sg, ug, wu, alt_label);
            }
            nb_cut++;
            if(verbose){
                fprintf(stderr, "[debug::%s] cut tip utg%.6d of utg%.6d\n", __func__, (int)(wu>>1)+1, (int)(vu>>1)+1);
            }
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d tips\n", __func__, nb_cut);
    }
    return nb_cut;
}

int hamt_ug_pop_simpleShortCut(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // FUNC
    /*
             --v1--
            /      \
        ...v0------v2---
        drop v0->v1 and v1->v2
    */
   // NOTE
   //    to this point, we expect few short tips,
   //    so this function doesn't tolerate tips.
    int verbose = 0;

    uint32_t v0, v1, v2, v[2], w, nv;
    asg_arc_t *av;
    asg_t *auxsg = ug->g;
    int nb_cut = 0;

    for (v0=0; v0<auxsg->n_seq*2; v0++){
        if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0, base_label)!=2){continue;}
        nv = asg_arc_n(auxsg, v0);
        av = asg_arc_a(auxsg, v0);
        int idx = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            v[idx++] = av[i].v;
            assert(idx<=2);
        }
        // check topo
        if (hamt_asgarc_util_countSuc(auxsg, v[0], 0, 0, base_label)!=1 && hamt_asgarc_util_countSuc(auxsg, v[1], 0, 0, base_label)!=1){
            continue;
        }
        if (hamt_asgarc_util_countPre(auxsg, v[0], 0, 0, base_label)==1 && hamt_asgarc_util_countPre(auxsg, v[1], 0, 0, base_label)==2){
            v1 = v[0];
            v2 = v[1];
        }else if (hamt_asgarc_util_countPre(auxsg, v[0], 0, 0, base_label)==2 && hamt_asgarc_util_countPre(auxsg, v[1], 0, 0, base_label)==1){
            v1 = v[1];
            v2 = v[0];
        }else{
            continue;
        }
        if (hamt_asgarc_util_countSuc(auxsg, v1, 0, 0, base_label)!=1){continue;}
        // check topo (the targeting on v2)
        if (hamt_asgarc_util_get_the_one_target(auxsg, v1, &w, 0, 0, base_label)<0){fprintf(stderr, "ERROR %s\n", __func__); fflush(stderr); exit(1);}
        if (w==v2){
            // check length: otherwise might break stuff
            if (ug->u.a[v1>>1].len>50000){continue;}
            // cut
            if (is_hard_drop){
                hamt_ug_arc_del(sg, ug, v0, v1, 1);
                hamt_ug_arc_del(sg, ug, v1, v2, 1);
            }else{
                hamt_ug_arc_del(sg, ug, v0, v1, 1);
                hamt_ug_arc_del(sg, ug, v1, v2, 1);
                hamt_ug_utg_softdel(sg, ug, v1, alt_label);
            }
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] cut route: utg%.6d utg%.6d utg%.6d\n", __func__, (int)(v0>>1)+1, (int)(v1>>1)+1, (int)(v2>>1)+1);
            }
        }
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d spots\n", __func__, nb_cut);
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return nb_cut;
}

int hamt_ug_oneutgCircleCut(asg_t *sg, ma_ug_t *ug, int base_label){
    // if a long unitig forms a circle on its own, and is only linked to other subgraphs
    // by 1 arc (the target also has only this one predecessor), cut it off regardless of coverage diff
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t vu, vu_, wu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=1){continue;}
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=2){continue;}
        if (hamt_asgarc_util_get_the_one_target(auxsg, vu, &vu_, 0, 0, base_label)<0){fprintf(stderr, "ERROR %s\n", __func__); fflush(stderr); exit(0);}
        if (vu_!=vu){continue;}  // check circling
        nv = asg_arc_n(auxsg, vu^1);
        av = asg_arc_a(auxsg, vu^1);
        int san = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if ((av[i].v^1)==vu){continue;}
            san = 1;
            wu = av[i].v;
        }
        assert(san);
        if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)!=1){continue;}  // backward branching
        if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)==0){  // spare if is long tip
            if (ug->u.a[wu>>1].len>100000){
                continue;
            }
        }  

        // topo check passed, cut
        hamt_ug_arc_del(sg, ug, vu^1, wu, 1);  // compile note: safe by assertion
        nb_cut++;
        if (verbose){
            fprintf(stderr, "[M::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
        }
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
    }
    return nb_cut;
}

int hamt_ug_oneutgCircleCut2(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //     1) if a long unitig forms a circle on its own, and is linked to other subgraphs
    //     by multiple arcs, check coverage diff. If the unitig's coverage is much higher/lower than ALL of
    //     the neighbouring nodes', cut it off.
    //     2) if the neighbours are one same unitig (i.e. circle), check coverage, if it's less than 2x diff,
    //        cut the inner linkage (i.e. merge the unitig in question and the neighbour unitig).
    // NOTE
    //     assume ug has its coverages collected.
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t vu, vu_, wu, nv, nv2;
    asg_arc_t *av, *av2;
    int nb_cut = 0;

    uint32_t idx;
    int cov0, cov_tmp, cov_min, cov_max, cov_type2=0;
    float ratio1, ratio2;


    for (vu=0; vu<auxsg->n_seq*2; vu++){
        cov_min=0xfffffff, cov_max=0; // reset
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label){continue;}
        if (ug->u.a[vu>>1].len<1000000){continue;}
        
        if (verbose>1){fprintf(stderr, "[debug::%s] at utg%.6d\n", __func__, (int)(vu>>1)+1);}
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<2){
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed suc\n", __func__);}
            continue;
        }
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (idx=0; idx<nv; idx++){
            if (av[idx].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[idx].v>>1]!=base_label) {continue;}
            if (av[idx].v==vu){break;}
        }
        if (idx==nv){
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed to target itself\n", __func__);}
            continue;
        }

        // the utg forms a circle. Check coverage
        cov0 = ug->utg_coverage[av[idx].v>>1];
        for (int i=0; i<nv; i++){
            if (av[i].del || (i==idx)) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            cov_tmp = ug->utg_coverage[av[i].v>>1];
            cov_min = cov_min>cov_tmp? cov_tmp : cov_min;
            cov_max = cov_max<cov_tmp? cov_tmp : cov_max;
        }

        if ( (cov0+15)<cov_min || (cov0*2)<cov_min || (cov_max*2)<cov0 || (cov_max+15)<cov0){
            // drop 1st side
            for (int i=0; i<nv; i++){
                if (av[i].del || (i==idx)) {continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
            }
            // drop 2nd side
            nv = asg_arc_n(auxsg, vu^1);
            av = asg_arc_a(auxsg, vu^1);
            for (int i=0; i<nv; i++){
                if (av[i].del || ((av[i].v^1)==vu) ){continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                hamt_ug_arc_del(sg, ug, vu^1, av[i].v, 1);
            }
            if (verbose){
                fprintf(stderr, "[debug::%s] (type1) dropped arcs for: utg%.6d\n", __func__, (int)(vu>>1)+1);
                fprintf(stderr, "[debug::%s]          cov %d, cov max %d, cov min%d\n", __func__, cov0, cov_max, cov_min);
            }
            nb_cut++;
        }else{
            if (verbose>1){fprintf(stderr, "[debug::%s]     failed cov diff\n", __func__);}
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
    }
    return nb_cut;
}


int hamt_ug_checknpop_oneMultiLeafSoapBubble(asg_t *sg, ma_ug_t *ug, uint32_t start0, 
                                             int base_label, int alt_label, int is_hard_drop)
{
    // FUNC
    //     (the naming isn't too expressive)
    //     check and pop the following topo:
    /*
                         v2
                         |
                        v2^1-----
                        /        \
                  v1~v1^1     v3~v3^1
                 /                |
    ...---v0~v0^1------------v4~v4^1-------v5~v5^1--...
        drop v0->v1 and v3->v4 (both directions)
        (note that only one end of vertex v2 is used. This structure can't be treated by complex bubble popping.)
    */
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu[6], wu[2], vu_tmp, nv;
    int san[2], suc;
    asg_arc_t *av;
    int idx = 0;

    // get v0, v1, v4
    if (hamt_asgarc_util_countSuc(auxsg, start0, 0, 0, base_label)!=2) {return 0;}
    nv = asg_arc_n(auxsg, start0);
    av = asg_arc_a(auxsg, start0);
    san[0] = 0; san[1] = 0;
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del){continue;}
        if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
        suc = hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label);
        if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)!=1 || suc>2 || suc==0){
            return 0;
        }
        wu[suc-1] = av[i].v;
        san[suc-1]++;
    }
    if (san[0]!=1 || san[1]!=1){return 0;}
    vu[0] = start0; // v0
    vu[1] = wu[0];  // v1
    vu[4] = wu[1];  // v4
    assert(hamt_asgarc_util_get_the_one_target(auxsg, vu[1], &vu[2], 0, 0, base_label)>=0);
    assert(hamt_asgarc_util_get_the_one_target(auxsg, vu[4], &vu[5], 0, 0, base_label)>=0);

    // check v5's topo
    if (hamt_asgarc_util_countPre(auxsg, vu[5], 0, 0, base_label)!=1){
        return 0;
    }
    
    // check v2's topo
    if (hamt_asgarc_util_countSuc(auxsg, vu[2], 0, 0, base_label)!=0){
        return 0;
    }
    if (hamt_asgarc_util_countPre(auxsg, vu[2], 0, 0, base_label)!=2){
        return 0;
    }

    // get v3
    nv = asg_arc_n(auxsg, vu[2]^1);
    av = asg_arc_a(auxsg, vu[2]^1);
    {
        uint32_t i;
        for (i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (av[i].v==(vu[1]^1)) {continue;}
            vu[3] = av[i].v;
            break;
        }
        assert(i!=nv);
    }
    if (hamt_asgarc_util_countPre(auxsg, vu[3], 0, 0, base_label)!=1 || hamt_asgarc_util_countSuc(auxsg, vu[3], 0, 0, base_label)!=1){
        return 0;
    }
    assert(hamt_asgarc_util_get_the_one_target(auxsg, vu[3], &vu_tmp, 0, 0, base_label)>=0);
    if ((vu_tmp^1)!=vu[4]){return 0;}

    // if we reach here, topo check is all clear. Drop the arcs
    if (is_hard_drop){
        hamt_ug_arc_del(sg, ug, vu[0], vu[1], 1);
        hamt_ug_arc_del(sg, ug, vu[3], vu[4]^1, 1);
    }else{
        hamt_ug_arc_del(sg, ug, vu[0], vu[1], 1);
        hamt_ug_arc_del(sg, ug, vu[3], vu[4]^1, 1);
        hamt_ug_utg_softdel(sg, ug, vu[1], alt_label);
        hamt_ug_utg_softdel(sg, ug, vu[2], alt_label);
        hamt_ug_utg_softdel(sg, ug, vu[3], alt_label);
    }
    if (verbose){fprintf(stderr, "[debug::%s] popped: start utg%.6d leaf utg%.6d\n", __func__, (int)(vu[0]>>1)+1, (int)(vu[2]>>1)+1);}
    return 1;
}


int hamt_ug_drop_midsizeTips(asg_t *sg, ma_ug_t *ug, int fold, int base_label){
    // FUNC
    //     treat the follow tip
    /*
               vu
              /
    --------wu----uu--------
        where either wu or uu is significantly longer than vu.
        Will cut the vu->wu arc. Won't send vu to the alternative collection.
    */
    asg_t *auxsg = ug->g;
    uint32_t wu, uu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    int r = fold<0? 5 : fold;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        // tip topo
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=1){continue;}

        // neighbour 1
        if (hamt_asgarc_util_get_the_one_target(auxsg, vu, &wu, 0, 0, base_label)<0){
            fprintf(stderr, "[E::%s] can't get target\n", __func__);
            continue;
        }
        if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)==0 || hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)!=2){continue;}

        // neighbour 2
        nv = asg_arc_n(auxsg, wu^1);
        av = asg_arc_a(auxsg, wu^1);
        uint32_t i;
        for (i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if ((av[i].v>>1)==(vu>>1)){continue;}
            uu = av[i].v;
            break;
        }
        if (i==nv){  // didn't found uu for some reason
            fprintf(stderr, "[E::%s] didn't found uu when should've\n", __func__);
            continue;
        }
        if (hamt_asgarc_util_countSuc(auxsg, uu, 0, 0, base_label)==0 || hamt_asgarc_util_countPre(auxsg, uu, 0, 0, base_label)!=1){continue;}  // compile note: is safe

        // check length
        if (ug->u.a[wu>>1].len>(ug->u.a[vu>>1].len*r) || ug->u.a[uu>>1].len>(ug->u.a[vu>>1].len*r)){
            // cut
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            hamt_ug_arc_del(sg, ug, wu^1, uu, 1);
            nb_cut++;
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d mid length tips.\n", __func__, nb_cut);
    }
    return nb_cut;
}

void hamt_ug_prectgTopoClean(asg_t *sg, 
                            const ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex,
                            int base_label, int alt_label, int is_hard_drop){
    // note: will redo the unitig graph
    // assume graph is clean
    // pop all simple bubbles, apparent bubbles, pentagon bubbles (popping doesn't depend on coverage)
    double startTime = Get_T();
    int verbose = 0;

    ma_ug_t *ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, base_label);
    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu, nv;
    asg_arc_t *av;
    int nb_pop = 0, round = 0;

    for (round=0; round<5; round++){
        if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, "prectg-TOPO", round);}

        nb_pop = 0;

        // inverted simple bubbles
        nb_pop += hamt_ug_pop_simpleInvertBubble(sg, ug, base_label, alt_label, is_hard_drop);

        // penta simple bubbles
        uint32_t penta_buf[3];
        int nb_penta_pop = 0;
        for (vu=0; vu<auxsg->n_seq*2; vu++){
            if (hamt_asgarc_util_checkSimplePentaBubble(auxsg, vu, penta_buf, base_label)>0){
                if (verbose){
                    fprintf(stderr, "[debug::%s] penta bubble pop: from utg%.6d , route utg%.6d utg%.6d utg%.6d\n", __func__,
                            (int)(vu>>1)+1, (int)(penta_buf[0]>>1)+1, (int)(penta_buf[1]>>1)+1, (int)(penta_buf[2]>>1)+1);
                }
                if (is_hard_drop){
                    hamt_ug_arc_del(sg, ug, penta_buf[0], penta_buf[1], 1);
                    hamt_ug_arc_del(sg, ug, penta_buf[1], penta_buf[2], 1);
                }else{
                    hamt_ug_arc_del(sg, ug, penta_buf[0], penta_buf[1], 1);
                    hamt_ug_arc_del(sg, ug, penta_buf[1], penta_buf[2], 1);
                    hamt_ug_utg_softdel(sg, ug, penta_buf[1], alt_label);
                }
                nb_pop++;
                nb_penta_pop++;
            }
        }
        if (verbose){
            fprintf(stderr, "[M::%s] popped %d penta simple bubbles (round %d)\n", __func__, nb_penta_pop, round);
        }

        // bi-link bubble chain
        int nb_bi_pop = 0;
        for (vu=0; vu<auxsg->n_seq*2; vu++){
            if (verbose){
                fprintf(stderr, "[debug::%s] bi-bubble, at utg%.6d\n", __func__, (int)(vu>>1)+1);
            }
            hamt_ug_util_popSimpleBiBubbleChain(sg, ug, vu, base_label);
            nb_bi_pop++;
            nb_pop++;
        }
        if (verbose){
            fprintf(stderr, "[M::%s] popped %d bi bubbles chains(round %d)\n", __func__, nb_bi_pop, round);
        }


        if (nb_pop==0){
            hamt_ug_destroy(ug);
            if (VERBOSE){
                fprintf(stderr, "[M::%s] (early termination at round %d)\n", __func__, round);
            }
            return;
        }else{
            asg_cleanup(sg);
            hamt_ug_regen(sg, &ug, coverage_cut, sources, ruIndex, base_label);
            auxsg = ug->g;
            if (VERBOSE){
                fprintf(stderr, "[M::%s] popped total %d locations (round %d).\n", __func__, nb_pop, round);
                fprintf(stderr, "[T::%s] took %0.2fs.\n\n", __func__, Get_T()-startTime);

            }
        }
    }
    hamt_ug_destroy(ug);

}

void hamt_ug_prectg_rescueShortCircuit(asg_t *sg, 
                                        ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, R_to_U* ruIndex,
                                        const ma_sub_t* coverage_cut, int base_label){
    // exp
    // addresses observations in the sheep datasets, where 2 apparent circles
    //   are joined by one relatively short unitig (not necessarily 1-read)
    //   if the ends of the circle actually overlap with each other, then ditch the short unitig and 
    //   recover that arc (which was dropped during transitive reduction)
    // FUNC
    //    recover some arc
    // NOTE
    //    use after basic topo clean has been done
    //        this function will only check very simple circles
    //    will regenerate the unitig graph, maybe several times
    //    hopefully this doesn't mess up with other stuffs?
    int verbose = 0;
    int nb_modified = 0;
    uint32_t start, end, start_v, end_v;  // start and end are meant to be unitig IDs (with dir); start_v and end_v are vertex IDs
    int l1, l2;
    uint32_t nv, wu, uu;
    asg_arc_t *av;
    int is_circular;

    for (int round=0; round<3; round++){
        if (verbose){
            fprintf(stderr, "[debug::%s] entered round %d\n", __func__, round);
        }
        ma_ug_t *ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, base_label);
        asg_t *auxsg = ug->g;

        // debug
        if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, "prectg-resShortCircuit", round);}

        for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
            if (verbose){
                fprintf(stderr, "[debug::%s] at utg%.6d\n", __func__, (int)(vu>>1)+1);
            }

            // collect start and end unitigs of the single path
            int length_bp = 0;
            l1=hamt_ugasg_util_findendSinglePath_allowSimpleBubble(ug, vu, &end, &is_circular, &length_bp, base_label);
            if (is_circular){
                if (verbose){
                    fprintf(stderr, "[debug::%s]     1st extension was circular\n", __func__);
                }
                continue;
            }
            l2=hamt_ugasg_util_findendSinglePath_allowSimpleBubble(ug, vu^1, &start, &is_circular, &length_bp, base_label);
            if (is_circular){
                if (verbose){
                    fprintf(stderr, "[debug::%s]     2nd extension was circular\n", __func__);
                }
                continue;
            }
            if (l1==0 && l2==0){
                length_bp = ug->u.a[vu>>1].len;
            }
            if (length_bp<1000000){
                if (verbose){
                    fprintf(stderr, "[debug::%s]     single pass too short (%d).\n", __func__, length_bp);
                }
                continue;
            }
            start^=1;
            if (hamt_asgarc_util_countSuc(auxsg, end, 0, 0, base_label)==0){continue;}
            if (hamt_asgarc_util_countPre(auxsg, start, 0, 0, base_label)==0){continue;}
            if (hamt_asgarc_util_countSuc(auxsg, end, 0, 0, base_label)!=1){
                if (l1>0){
                    fprintf(stderr, "something wrong 1a%s, vu is utg%.6d, start utg%.6d, end utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(start>>1)+1, (int)(end>>1)+1);
                }
                continue;
            }
            if (hamt_asgarc_util_countPre(auxsg, start, 0, 0, base_label)!=1){
                if (l2>0){
                    fprintf(stderr, "something wrong 1b%s, vu is utg%.6d, start utg%.6d, end utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(start>>1)+1, (int)(end>>1)+1);
                }
                continue;
            }

            // check if after stepping 1 step, the end unitg will find the start unitig
            if (verbose){
                fprintf(stderr, "[debug::%s]     got to check if there was any overlap\n", __func__);
            }
            if (hamt_asgarc_util_get_the_one_target(auxsg, end, &wu, 0, 0, base_label)<0){
                fprintf(stderr, "something wrong 2%s\n", __func__);
                continue;

            }
            if (ug->u.a[wu>>1].len>100000){  // require the utg in question to be short
                if (verbose){
                    fprintf(stderr, "[debug::%s]     mid utg too long\n", __func__);
                }
                continue;
            }
            nv = asg_arc_n(auxsg, wu);
            av = asg_arc_a(auxsg, wu);
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del){continue;}
                if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                uu = av[i].v;
                if (uu==start){  // check if there WAS an arc conneting end and start
                    // get vertices
                    if (end&1){
                        end_v = (ug->u.a[end>>1].start)^1;
                    }else{
                        end_v = (ug->u.a[end>>1].end)^1;
                    }
                    if (start&1){
                        start_v = ug->u.a[start>>1].end;
                    }else{
                        start_v = ug->u.a[start>>1].start;
                    }
                    if (verbose){
                        fprintf(stderr, "[debug::%s]     end vertex is %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, end_v>>1), Get_NAME(R_INF, end_v>>1));
                        fprintf(stderr, "[debug::%s]     start vertex is %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, start_v>>1), Get_NAME(R_INF, start_v>>1));
                    }
                    // check if overlap ever existed
                    if (hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, 5)>0){
                        // drop old links
                        hamt_ug_arc_del(sg, ug, vu, wu, 1);
                        hamt_ug_arc_del(sg, ug, wu, uu, 1);
                        // log
                        nb_modified++;
                        if (verbose){
                            fprintf(stderr, "[debug::%s] treated utg%.6d, utg%.6d, utg %.6d \n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
                        }
                    }
                    break;  // since there's no multi arc, we only expect uu==start once.
                }
            }
        }
        // hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, 210+round);
        hamt_ug_destroy(ug);
        asg_cleanup(sg);
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] modified %d spots\n", __func__, nb_modified);
    }
}

void hamt_ug_prectg_rescueShortCircuit_simpleAggressive(asg_t *sg, ma_ug_t *ug, 
                                                        ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,
                                                        const ma_sub_t* coverage_cut, 
                                                        int base_label){
    // FUNC
    //     simpler verssion of hamt_ug_prectg_rescueShortCircuit
    //     if a subgraph contains v0->v1->v0 loop (v1 can have other branches but v0 dont),
    //       and v1 is 1) much shorter than v0, and 2) has a different coverage than v0
    //     then check if there was any overlap between the start and end of v0. If so,
    //       cut the arc between v0 and v1, then recover the v0->v0 arc.
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t nv, wu, uu, v, w;
    asg_arc_t *av;
    int is_passed=0, nb_treat = 0;
    float diff = 0;
    ma_hit_t *h;

    for (uint32_t v1=0; v1<auxsg->n_seq*2; v1++){
        if (base_label>=0 && auxsg->seq_vis[v1>>1]!=base_label){continue;}
        if (hamt_asgarc_util_countPre(auxsg, v1, 0, 0, base_label)==0 || hamt_asgarc_util_countSuc(auxsg, v1, 0, 0, base_label)==0){continue;}
        nv = asg_arc_n(auxsg, v1);
        av = asg_arc_a(auxsg, v1);

        // check v1's targets
        for (uint32_t i=0; i<nv; i++){
            is_passed = 0;

            if (av[i].del) {continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            
            // v0 shall target v1 to qualify
            wu = av[i].v;
            if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)==1 && 
                hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)==1 &&
                hamt_asgarc_util_get_the_one_target(auxsg, wu, &uu, 0, 0, base_label)>=0){
                if (uu==v1){
                    // check coverage (don't touch if v0 and v1 have similar coverage)
                    diff = (float)ug->utg_coverage[v1>>1]/ug->utg_coverage[wu>>1];
                    if (diff>=0.33 && diff<=3){
                        continue;
                    }
                    // check if there's ever an overlap 
                    v = ug->u.a[wu>>1].start;
                    w = ug->u.a[wu>>1].end^1;
                    if (hamt_ug_recover_ovlp_if_existed(sg, ug, w, v, sources, coverage_cut, 5)>0){
                        if (verbose){
                            fprintf(stderr, "[debug::%s] treated utg%.6d\n", __func__, (int)(wu>>1)+1);
                        }
                        // cut old arcs
                        hamt_ug_arc_del(sg, ug, v1, wu, 1);
                        nb_treat++;
                    }else{
                        if (verbose){
                            fprintf(stderr, "[debug::%s] utg%.6d passed topo, but had no overlap\n", __func__, (int)(wu>>1)+1);
                        }
                    }
                }
            }
        }
    }

    if (nb_treat){
        // hamt_ug_cleanup_arc_by_labels(sg, ug);
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_treat);
    }

}

void hamt_ug_prectg_rescueLongUtg(asg_t *sg, 
                                    ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,R_to_U* ruIndex,
                                    const ma_sub_t* coverage_cut)
{
    // experimental
    //    if a long unitig is not connected and non-circular,
    //    try the ends and add a link if there was an overlap in the paf
    //    Obviously more aggressive than hamt_ug_prectg_rescueShortCircuit, which required more topo checks. 
    int verbose = 1;

    ma_ug_t *ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, 0);
    asg_t *auxsg = ug->g;
    uint32_t start_v, end_v;
    int nb_treated = 0;
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq, 1);  // easier than fixing auxsg

    if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, "prectg-rescueLongUtg",0);}

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){  
        if (color[vu>>1]){continue;}
        if (ug->u.a[vu>>1].len<1000000){
            continue;
        }
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, -1)!=0 || hamt_asgarc_util_countPre(auxsg, vu, 0, 0, -1)!=0){
            continue;
        }

        // check and add link if found
        if (verbose) {fprintf(stderr, "rescueLong\tat utg%.6d\n", (int)(vu>>1)+1);}
        start_v = ug->u.a[vu>>1].start;
        end_v = ug->u.a[vu>>1].end^1;
        if (hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, 20)>0){
            nb_treated++;
            color[vu>>1] = 1;
            if (verbose){
                fprintf(stderr, "[debug::%s] added circle link for utg%.6d \n", __func__, (int)(vu>>1)+1);
            }
        }        
    }

    // hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, 221);
    if (nb_treated){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
    }
    hamt_ug_destroy(ug);

    free(color);
    if (VERBOSE){
        fprintf(stderr, "[M::%s] added %d links.\n", __func__, nb_treated);
    }
}


int hamt_ug_prectg_resolve_complex_bubble(asg_t *sg, ma_ug_t *ug, 
                                          int base_label, int alt_label, int is_hard_drop){
    int verbose = 2;

    // ma_ug_t *ug = ma_ug_gen(sg);
    asg_t *auxsg = ug->g;
    uint32_t vu, nv, end, wu;
    asg_arc_t *av;
    int nb_cut = 0, nb_tangles = 0;

    vecu64_t vec, vec_seen;
    vecu64_init(&vec);
    vecu64_init(&vec_seen);

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        // (expects vu to be the start of a complex bubble)
        if (hamt_ug_check_complexBubble(sg, ug, 500, vu, &end, base_label)>0){
            vecu64_push(&vec, ( ((uint64_t)vu)<<32 )|( (uint64_t)end ) );
            nb_tangles++;
        }
    }

    uint64_t key, key_rev;
    if (nb_tangles){
        for (int i=0; i<vec.n; i++){  // TODO: this is ugly
            key = vec.a[i];
            vu = (uint32_t) (key>>32);
            end = (uint32_t)key;
            key_rev = ( ((uint64_t)end^1)<<32 )|( (uint64_t)(vu^1) );

            if (vecu64_is_in_vec(&vec_seen, key_rev)){
                if (verbose){fprintf(stderr, "[debug::%s] skip utg%.6d -> utg%.6d bc seen\n", __func__, (int)(vu>>1)+1, (int)(end>>1)+1);}
                continue;
            }
            
            if (!vecu64_is_in_vec(&vec, key_rev)){
                if (verbose){
                    fprintf(stderr, "[W::%s] only identified in one direction for utg%.6d and utg%.6d, will not treat\n", 
                                    __func__, (int)(vu>>1)+1, (int)(end>>1)+1);
                }
                continue;
            }
            if (verbose){
                fprintf(stderr, "[W::%s] treatment for utg%.6d and utg%.6d\n", 
                                    __func__, (int)(vu>>1)+1, (int)(end>>1)+1);
            }
            nb_cut += hamt_ug_pop_complexBubble(sg, ug, vu, end, base_label, alt_label, is_hard_drop);   
            vecu64_push(&vec_seen, key);
        }
    }
    vecu64_destroy(&vec);
    vecu64_destroy(&vec_seen);

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        hamt_ug_cleanup_arc_by_labels(sg, ug);
    }

    if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, "resolveTangle_middle", 0);}

    // cut tips
    int nb_cut_tip = 1, nb_cut_tiny_circle = 1;
    int tip_round = 0;
    while (nb_cut_tip>0 || nb_cut_tiny_circle>0){
        nb_cut_tip = hamt_ugasg_cut_shortTips(sg, ug, base_label, alt_label, is_hard_drop);
        nb_cut_tip += hamt_ug_cut_shortTips_arbitrary(sg, ug, 30000, base_label);
        nb_cut_tiny_circle = hamt_ug_pop_tinyUnevenCircle(sg, ug, base_label, alt_label, is_hard_drop);
        if (VERBOSE){
            fprintf(stderr, "[M::%s] topo clean up, dropped %d tips, %d tiny uneven circles (round %d)\n", __func__, nb_cut_tip, nb_cut_tiny_circle, tip_round);
        }
        tip_round++;
        hamt_ug_cleanup_arc_by_labels(sg, ug);
        if (tip_round>10){break;}  // waterproof just in case
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d arcs for %d complex bubbles\n", __func__, nb_cut, nb_tangles);
    }
    return nb_cut;
}

int hamt_ug_resolve_oneMultiLeafSoapBubble(asg_t *sg, ma_ug_t *ug, 
                                           int base_label, int alt_label, int is_hard_drop){
    int nb_cut = 0;
    for (uint32_t vu=0; vu<ug->g->n_seq*2; vu++){
        nb_cut += hamt_ug_checknpop_oneMultiLeafSoapBubble(sg, ug, vu, base_label, alt_label, is_hard_drop);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
    }
    return nb_cut;
}

int hamt_ug_resolve_small_multileaf_with_covcut(asg_t *sg, ma_ug_t *ug, int max_length, int fold, int base_label){
    // FUNC
    //     more of a post-processing multileaf cleaning
    int nb_cut = 0;
    asg_t *auxsg = ug->g;
    uint32_t wu, nv;
    asg_arc_t *av;

    int ml = max_length<0? 50000 : max_length;
    int r = fold<0? 3 : fold;
    int handle_cov;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        if (ug->u.a[vu>>1].len>max_length){continue;}
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<2){continue;}
        handle_cov = ug->utg_coverage[vu>>1]==0? 1 : ug->utg_coverage[vu>>1];

        // drop arcs with coverage diff
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            wu = av[i].v;
            if (ug->utg_coverage[wu>>1]>=(r*handle_cov) || (r*handle_cov)>=ug->utg_coverage[wu>>1]){
                // cut
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                nb_cut++;
            }
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d arcs\n", __func__, nb_cut);
    }
    return nb_cut;
}


int hamt_ug_treatBifurcation_hapCovCut(asg_t *sg, ma_ug_t *ug, float covdiff_ratio, float haplo_ratio, 
                                        ma_hit_t_alloc *reverse_sources,
                                        int base_label, int alt_label){
    // FUNC
    //    if vertex vu has 2 targets and they appeart to form a pair of haplotigs,
    //    check coverage, if one of the targets appears to have a compatible coverage with vu,
    //    cut the other arc.
    // RETURN
    //    nb_cut
    int verbose = 0;

    asg_t *auxsg = ug->g;
    int nb_cut = 0, san, idx;
    uint32_t wu[2], nv;
    asg_arc_t *av;

    float r1, r2;  // for checking coverage diff

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        // specify topo: forward bifurcation
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=2){continue;}
        san = 0;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            if (san>=2){
                fprintf(stderr, "[E::%s] sancheck failed (1); skip and continue. (san: %d)\n", __func__, san);
            }
            wu[san] = av[i].v;
            san++;
        }
        if (san!=2){
            fprintf(stderr, "[E::%s] sancheck failed (2); skip and continue. (san: %d)\n", __func__, san);
            continue;
        }

        // only treat long unitigs
        if (ug->u.a[vu>>1].len<100000 || ug->u.a[wu[0]>>1].len<100000 || ug->u.a[wu[1]>>1].len<100000){continue;}


        if (verbose) {fprintf(stderr, "[debug::%s] at utg%.6d, targeting utg%.6d and utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);}

        // check coverage
        // (don't check haplo if coverage won't fulfill the criteira anyway)
        if (ug->utg_coverage[vu>>1] < ug->utg_coverage[wu[0]>>1]){
            r1 = (float)ug->utg_coverage[vu>>1]/ug->utg_coverage[wu[0]>>1];
        }else{
            r1 = (float)ug->utg_coverage[wu[0]>>1]/ug->utg_coverage[vu>>1];
        }
        if (ug->utg_coverage[vu>>1] < ug->utg_coverage[wu[1]>>1]){
            r2 = (float)ug->utg_coverage[vu>>1]/ug->utg_coverage[wu[1]>>1];
        }else{
            r2 = (float)ug->utg_coverage[wu[1]>>1]/ug->utg_coverage[vu>>1];
        }
        if (r1<covdiff_ratio && r2<covdiff_ratio){
            if (verbose) {fprintf(stderr, "[debug::%s]     didn't pass coverage check\n", __func__);}
            continue;
        }
        if (r1>=covdiff_ratio && r2>=covdiff_ratio){
            if (verbose) {fprintf(stderr, "[debug::%s]     both passed coverage check, don't touch to be safe\n", __func__);}
            continue;
        }
        idx = r1>r2? 1:0;  // wu[idx] is the target with more coverage diff
        if (verbose) {fprintf(stderr, "[debug::%s]     coverage check: might remove utg%.6d\n", __func__, (int)(wu[idx]>>1)+1);}
        
        // check if the two form a pair of haplotigs
        if (hamt_check_diploid(ug, wu[0], wu[1], haplo_ratio, reverse_sources)>0){
            // cut
            hamt_ug_arc_del(sg, ug, vu, wu[idx], 1);
            // hamt_ug_utg_softdel(sg, ug, wu[idx], alt_label);  // DON'T!
            nb_cut++;
            if (verbose) {fprintf(stderr, "[debug::%s]     CUT\n", __func__);}
        }else{
            if (verbose) {fprintf(stderr, "[debug::%s]     didn't pass hap check\n", __func__);}
        }
        
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
    }
    return nb_cut;
}




int hamt_ug_basic_topoclean(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    int nb_cut = 0;
    nb_cut += hamt_ug_pop_bubble(sg, ug, base_label, alt_label, is_hard_drop);  // note: include small tip cutting
    // hamt_asgarc_ugTreatMultiLeaf(sg, ug, 50000);  // note: doesn't protect obvious end-of-path tips
    nb_cut += hamt_ug_pop_miscbubble(sg, ug, base_label);
    nb_cut += hamt_ug_pop_simpleInvertBubble(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_pop_miscbubble_aggressive(sg, ug, base_label);

    nb_cut += hamt_ug_pop_terminalSmallTip(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_pop_tinyUnevenCircle(sg, ug, base_label, alt_label, is_hard_drop);

    nb_cut += hamt_ug_pop_simpleShortCut(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_oneutgCircleCut(sg, ug, base_label);
    return nb_cut;
}

int hamt_ug_basic_topoclean_simple(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    int nb_cut = 0;
    nb_cut += hamt_ug_pop_bubble(sg, ug, base_label, alt_label, is_hard_drop);  // note: include small tip cutting
    nb_cut += hamt_ug_pop_terminalSmallTip(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_pop_tinyUnevenCircle(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_pop_simpleShortCut(sg, ug, base_label, alt_label, is_hard_drop);
    nb_cut += hamt_ug_oneutgCircleCut(sg, ug, base_label);
    return nb_cut;
}

#if 0
//
//  NOTE
//     Form chimeric haplotype for low coverage regions that would otherwise end up with 2 or more unitigs.
//     Has performance issues (could be resolve by using kt_pipeline or rewrite),
//       and seems to matter very little for most cases: if a speicies has too many this 
//       kind of spots, the whole assembly will still be fragmented even if we rescue - too many real gaps.
//     TODO - come back to this if more examples emerge. (Dec 9 2020; example: zymo(std), Candida)
//
int hamt_ug_rescueLowCovHapGap_simple(asg_t *sg, ma_ug_t *ug, 
                              ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, 
                              const ma_sub_t* coverage_cut, uint64_t *readLen){
    // Join two low coverage unitigs if a het arc is available
    // WILL modify sources.
    int verbose = 1;

    int nb_treated = 0;
    asg_t *auxsg = ug->g;
    uint32_t start, end;
    int ret;
    ma_hit_t *handle, handle_rev;

    vecu64_t seen;
    vecu64_init(&seen);
    vecu32_t treated;
    vecu32_init(&treated);
    uint64_t key;

    // collect unitig end info
    // -1 if unitig is not a tip
    // 1 if the end is simple (won't have more than 1 target candidate even if consider all overlaps)
    // 0 otherwise
    int32_t *read2utg = (int32_t*)calloc(sg->n_seq, sizeof(int32_t));  // TODO: use ruIndex?
    memset(read2utg, -1, sg->n_seq*sizeof(int32_t));
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) {continue;}
        for (int i=0; i<ug->u.a[vu>>1].n; i++){
            read2utg[ug->u.a[vu>>1].a[i]>>33] = vu>>1;
        }
    }
    int8_t *buf = (int8_t*)calloc(auxsg->n_seq, 1);
    memset(buf, 0, auxsg->n_seq*1);
    uint32_t prv_utg_ID;
    int breakout=0;
    uint32_t v;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        breakout = 0;
        if (vu&1) {continue;}
        fprintf(stderr, "marksimple\tutg%.6d\t", (int)(vu>>1)+1);

        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, -1)>0){
            buf[vu>>1] = -1; 
            fprintf(stderr, "-1\n");
            continue;
        }
        v = ug->u.a[vu>>1].end^1;
        // vu has no target in the current unitig graph; check overlaps
        prv_utg_ID = vu>>1;
        for (int i=0; i<sources[v>>1].length; i++){
            if (read2utg[sources[v>>1].buffer[i].tn]==-1){continue;}
            if (prv_utg_ID==(vu>>1) && read2utg[sources[v>>1].buffer[i].tn]!=(vu>>1)) {
                prv_utg_ID = read2utg[sources[v>>1].buffer[i].tn];
                fprintf(stderr, " (set: %.*s)\n", (int)Get_NAME_LENGTH(R_INF, sources[v>>1].buffer[i].tn), 
                                                             Get_NAME(R_INF, sources[v>>1].buffer[i].tn));
                continue;
            }
            if (read2utg[sources[v>>1].buffer[i].tn]!=prv_utg_ID){
                buf[vu>>1] = 0;
                fprintf(stderr, "0\n");
                fprintf(stderr, " (violation: %.*s)\n", (int)Get_NAME_LENGTH(R_INF, sources[v>>1].buffer[i].tn), 
                                                             Get_NAME(R_INF, sources[v>>1].buffer[i].tn));
                breakout = 1;
                // break;
            }
        }
        if (!breakout) {
            // check het ovlp
            for (int i=0; i<reverse_sources[v>>1].length; i++){
                if (read2utg[reverse_sources[v>>1].buffer[i].tn]==-1){continue;}
                if (prv_utg_ID==(vu>>1) && read2utg[reverse_sources[v>>1].buffer[i].tn]!=(vu>>1)) {
                    prv_utg_ID = read2utg[reverse_sources[v>>1].buffer[i].tn];
                    fprintf(stderr, " (set: %.*s)\n", (int)Get_NAME_LENGTH(R_INF, reverse_sources[v>>1].buffer[i].tn), 
                                                             Get_NAME(R_INF, reverse_sources[v>>1].buffer[i].tn));
                    continue;
                }
                if (read2utg[reverse_sources[v>>1].buffer[i].tn]!=prv_utg_ID){
                    buf[vu>>1] = 0;
                    breakout = 1;
                    fprintf(stderr, "0\n");
                    fprintf(stderr, " (violation: %.*s)\n", (int)Get_NAME_LENGTH(R_INF, reverse_sources[v>>1].buffer[i].tn), 
                                                             Get_NAME(R_INF, reverse_sources[v>>1].buffer[i].tn));
                    // break;
                }
            }
            if (!breakout){
                buf[vu>>1] = 1;
                fprintf(stderr, "1\n");
            }
        }
    }
    free(read2utg);

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) {continue;}// (only check vu's targets here, the other direction will be addressed by the iteration)
        // ovlp topo check: 
        if (buf[vu>>1]!=1){continue;}  // unitig's end might target more than one other unitigs

        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)!=0){continue;}
        for (uint32_t wu=0; wu<auxsg->n_seq*2; wu++){
            if ((vu>>1)==(wu>>1)){continue;}
            key = vu<wu? (((uint64_t)vu>>1)<<32 | (wu>>1)) : (((uint64_t)wu>>1)<<32 | (vu>>1)) ;
            if (vecu64_is_in_vec(&seen, key)){continue;}
            vecu64_push(&seen, key);
            if (verbose) {fprintf(stderr, "[debug::%s] at utg%.6d vs utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);}

            // wu forward
            if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, 0)==0){
                start = ug->u.a[vu>>1].end^1;
                end = ug->u.a[wu>>1].start;
                if ( (hamt_ovlp_read_coverage_nbreads(sources, start>>1, 0, readLen[start>>1]) + 
                      hamt_ovlp_read_coverage_nbreads(reverse_sources, start>>1, 0, readLen[end>>1])) <= 10){
                    // if (verbose) {fprintf(stderr, "cov check passed\n");}
                    ret = hamt_ug_recover_ovlp_if_existed(sg, ug, start, end, reverse_sources, coverage_cut, 10);
                    if (ret>0) {
                        nb_treated++;
                        // modify sources
                        handle = get_specific_overlap_handle(reverse_sources, start>>1, end>>1);
                        assert(handle>0);
                        add_ma_hit_t_alloc(&sources[start>>1], handle);
                        set_reverse_overlap(&handle_rev, handle);
                        add_ma_hit_t_alloc(&sources[end>>1], &handle_rev);
                        if (verbose){
                            fprintf(stderr, "wu backward sucess %d, %.*s - %.*s\n", ret, 
                                                    (int)Get_NAME_LENGTH(R_INF, start>>1), Get_NAME(R_INF, start>>1),
                                                    (int)Get_NAME_LENGTH(R_INF, end>>1), Get_NAME(R_INF, end>>1));
                        }
                        vecu32_push(&treated, wu>>1);
                        continue;  
                    }
                }
            }

            // wu backward
            if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, 0)==0){
                start = ug->u.a[vu>>1].end^1;
                end = ug->u.a[wu>>1].end;
                if ( (hamt_ovlp_read_coverage_nbreads(sources, start>>1, 0, readLen[start>>1]) + 
                      hamt_ovlp_read_coverage_nbreads(reverse_sources, start>>1, 0, readLen[end>>1])) > 10) {continue;}  // not low coverage
                // if (verbose) {fprintf(stderr, "cov check passed\n");}
                ret = hamt_ug_recover_ovlp_if_existed(sg, ug, start, end, reverse_sources, coverage_cut, 10);
                if (ret>0) {
                    nb_treated++;
                    // modify sources
                    handle = get_specific_overlap_handle(reverse_sources, start>>1, end>>1);
                    assert(handle>0);
                    add_ma_hit_t_alloc(&sources[start>>1], handle);
                    set_reverse_overlap(&handle_rev, handle);
                    add_ma_hit_t_alloc(&sources[end>>1], &handle_rev);
                    if (verbose){
                        fprintf(stderr, "wu backward sucess %d, %.*s - %.*s\n", ret, 
                                                (int)Get_NAME_LENGTH(R_INF, start>>1), Get_NAME(R_INF, start>>1),
                                                (int)Get_NAME_LENGTH(R_INF, end>>1), Get_NAME(R_INF, end>>1));
                    }  
                    vecu32_push(&treated, wu>>1);
                }
                
            }
        }        
    }
    vecu64_destroy(&seen);
    vecu32_destroy(&treated);
    if (nb_treated){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
        asg_arc_del_multi(sg);

        asg_cut_tip(sg, asm_opt.max_short_tip);
    }
    if (VERBOSE){
        fprintf(stderr, "[debug::%s] treated %d spots\n", __func__, nb_treated);
    }
    return nb_treated;

}


int hamt_ughit_rescueLowCovHapGap_core(asg_t *sg, ma_ug_t *ug,
                                            ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, const ma_sub_t* coverage_cut,
                                            long long n_read, uint64_t *readLen,
                                            long i_v, long i_w, int read_cov_threshold){
    // RETRURN
    //     1 if the arc should be recovered
    //     0 otherwise
    // SIDE EFFECT
    //     will copy overlap and set .del status accordingly for source and reverse_sources
    // NOTE
    //     Currently when looking for bridge read, this function only pick the 1st
    //       qualified candidate and doesn't check if there's more. This shouldn't
    //       matter too much?
    int verbose = 1;
    asg_t *auxsg = ug->g;

    if (hamt_ovlp_read_coverage_nbreads(sources, i_v, 0, readLen[i_v])>read_cov_threshold) {return 0;}  // coverage not low enough
    if (hamt_ovlp_read_coverage_nbreads(sources, i_w, 0, readLen[i_w])>read_cov_threshold) {return 0;}  // coverage not low enough

    uint32_t tn, qn;
    ma_hit_t *handle1, *handle2;
    int found=0, is_reverse_1, is_reverse_2;
    int iov_1, iov_2;

    ma_hit_t *placeholder1, *placeholder2;

    if (verbose) {fprintf(stderr, "[debug::%s] %.*s vs %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, i_v), Get_NAME(R_INF, i_v),
                                                                          (int)Get_NAME_LENGTH(R_INF, i_w), Get_NAME(R_INF, i_w));}
    // sancheck: the two reads to be bridged shall not overlap with each other
    if (does_ovlp_ever_exist(sources, i_v<<1, (uint32_t)i_w<<1, &placeholder1, &placeholder2) ||
        does_ovlp_ever_exist(reverse_sources, i_v<<1, (uint32_t)i_w<<1, &placeholder1, &placeholder2) ){
            if (verbose) {fprintf(stderr, "[debug::%s]     failed, bc two reads overlapped with each other\n", __func__);}
            return 0;
    }

    // search if there's a read linking query read (v) and target read (w)
    // (note: one of the overlap could be intra-haplotype, since phasing requires a few (3) reads.
    //        Therefore here we check both inter+inter and intra+inter cases.)
    for (iov_1=0; iov_1<reverse_sources[i_v].length; iov_1++){
        handle1 = &reverse_sources[i_v].buffer[iov_1];
        if (handle1->del){continue;}
        tn = handle1->tn;
        qn = (uint32_t)(handle1->qns>>32);
        if (verbose) {fprintf(stderr, "[debug::%s]     bridge is %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn));}
        for (iov_2=0; iov_2<reverse_sources[i_w].length; iov_2++){
            handle2 = &reverse_sources[i_w].buffer[iov_2];
            if (handle2->del){continue;}
            if ( ((uint32_t)(handle2->qns>>32))!=qn && handle2->tn==tn){
                if (hamt_ovlp_read_coverage_nbreads(sources, tn, 0, readLen[tn]) + 
                    hamt_ovlp_read_coverage_nbreads(reverse_sources, tn, 0, readLen[tn]) > read_cov_threshold) {continue;}
                found = 1; is_reverse_1=1; is_reverse_2=1;
                break;
            }
        }
        if (found) {break;}
        for (iov_2=0; iov_2<sources[i_w].length; iov_2++){
            handle2 = &sources[i_w].buffer[iov_2];
            if (handle2->del){continue;}
            if (handle2->tn==tn){
                if (hamt_ovlp_read_coverage_nbreads(sources, tn, 0, readLen[tn]) + 
                    hamt_ovlp_read_coverage_nbreads(reverse_sources, tn, 0, readLen[tn]) > read_cov_threshold) {continue;}
                found = 1; is_reverse_1=1; is_reverse_2=0;
                break;
            }
        }
        if (found) {break;}
    }
    if (!found){
        for (iov_1=0; iov_1<sources[i_v].length; iov_1++){
            handle1 = &sources[i_v].buffer[iov_1];
            if (handle1->del){continue;}
            tn = handle1->tn;
            qn = (uint32_t) (handle1->qns>>32);
            for (iov_2=0; iov_2<reverse_sources[i_w].length; iov_2++){
                handle2 = &reverse_sources[i_w].buffer[iov_2];
                if (handle2->del){continue;}
                if ( ((uint32_t)(handle2->qns>>32))!=qn && handle2->tn==tn){
                    if (hamt_ovlp_read_coverage_nbreads(sources, tn, 0, readLen[tn]) + 
                        hamt_ovlp_read_coverage_nbreads(reverse_sources, tn, 0, readLen[tn]) > read_cov_threshold) {continue;}
                    found = 1; is_reverse_1=0; is_reverse_2=1;
                    break;
                }
            }
            if (found){break;}
        }
    }

    if (found){  // found the bridging (also low-cov) read
        ma_hit_t handle_tmp;

        // add the overlaps to intra-haplotype collection && mark het counterparts as deleted
        if (is_reverse_1) {
            add_ma_hit_t_alloc(&sources[i_v], handle1);
            set_reverse_overlap(&handle_tmp, handle1);
            add_ma_hit_t_alloc(&sources[tn], &handle_tmp);
            handle1->del = 1;
        }
        if (is_reverse_2) {
            add_ma_hit_t_alloc(&sources[i_w], handle2);
            set_reverse_overlap(&handle_tmp, handle2);
            add_ma_hit_t_alloc(&sources[tn], &handle_tmp);
            handle2->del = 1;
        }

        fprintf(stderr, "hapgap closes: %.*s and %.*s\n", (int)Get_NAME_LENGTH(R_INF, i_v), Get_NAME(R_INF, i_v),
                                                            (int)Get_NAME_LENGTH(R_INF, i_w), Get_NAME(R_INF, i_w));

        // TODO: maybe don't need sorting...
        ma_hit_sort_tn(sources[i_v].buffer, sources[i_v].length);
        ma_hit_sort_tn(sources[i_w].buffer, sources[i_w].length);

        // modify graph
        if (is_reverse_1) {hamt_ug_recover_ovlp_if_existed_core(sg, ug, i_v<<1, tn<<1, reverse_sources, coverage_cut, 1);}  // note: shift is dummy, the last bit doesn't matter
        else{hamt_ug_recover_ovlp_if_existed_core(sg, ug, i_v<<1, tn<<1, sources, coverage_cut, 1);}
        if (is_reverse_2){hamt_ug_recover_ovlp_if_existed_core(sg, ug, i_w<<1, tn<<1, reverse_sources, coverage_cut, 1);}
        else{hamt_ug_recover_ovlp_if_existed_core(sg, ug, i_w<<1, tn<<1, sources, coverage_cut, 1);}

        return 1;
    }else{
        return 0;
    }
}
void hamt_ughit_rescueLowCovHapGap(asg_t *sg, ma_ug_t *ug, 
                                         ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, const ma_sub_t* coverage_cut,
                                         long long n_read, uint64_t *readLen, int read_cov_threshold){
    // FUNC
    //     Sometimes a low coverage region need to use overlaps between different haplotypes to
    //       not break into multiple unitigs. Although joining them isn't ideal, it's still preferable.
    //     This function treats such spots, where the gap can be closed by one read.
    // NOTE
    //     Will move the inter-haplotype overlap(s) to the intra-haplotype set as needed.
    // TODO
    //     Log ^this info. Might be useful for users.
    double startTime = Get_T();
    int verbose = 1;

    int nb_treated = 0;
    asg_t *auxsg = ug->g;
    uint32_t vid, wid;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1){continue;}  // only need to check utg in one direction
        if (verbose) {fprintf(stderr, "[debug::%s] at utg%.6d\n", __func__, (int)(vu>>1)+1);}
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, -1)==0){
            vid = ug->u.a[vu>>1].start>>1;
            for (uint32_t wu=0; wu<auxsg->n_seq*2; wu++){
                // the two unitigs to be joined shall not be a pair of haplotigs
                if (hamt_check_diploid(ug, vu, wu, 0.7, sources) || hamt_check_diploid(ug, vu, wu, 0.7, reverse_sources)){continue;}
                if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, -1)==0){  // left side of wu
                    wid = ug->u.a[wu>>1].start>>1;
                    nb_treated += hamt_ughit_rescueLowCovHapGap_core(sg, ug, sources, reverse_sources, coverage_cut, 
                                                       n_read, readLen, vid, wid, read_cov_threshold);          
                }
                if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, -1)==0){  // right side of wu
                    wid = ug->u.a[wu>>1].end>>1;
                    nb_treated += hamt_ughit_rescueLowCovHapGap_core(sg, ug, sources, reverse_sources, coverage_cut, 
                                                       n_read, readLen, vid, wid, read_cov_threshold);
                }
            }
        }
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, -1)==0){
            vid = ug->u.a[vu>>1].end>>1;
            for (uint32_t wu=0; wu<auxsg->n_seq*2; wu++){
                // the two unitigs to be joined shall not be a pair of haplotigs
                if (hamt_check_diploid(ug, vu, wu, 0.7, sources) || hamt_check_diploid(ug, vu, wu, 0.7, reverse_sources)){continue;}
                if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, -1)==0){  // left side of wu
                    wid = ug->u.a[wu>>1].start>>1;
                    nb_treated += hamt_ughit_rescueLowCovHapGap_core(sg, ug, sources, reverse_sources, coverage_cut, 
                                                       n_read, readLen, vid, wid, read_cov_threshold);
                }
                if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, -1)==0){  // right side of wu
                    wid = ug->u.a[wu>>1].end>>1;
                    nb_treated += hamt_ughit_rescueLowCovHapGap_core(sg, ug, sources, reverse_sources, coverage_cut, 
                                                       n_read, readLen, vid, wid, read_cov_threshold);
                }
            }
        }

    }

    if (nb_treated){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] done, %d spots involved. Took %0.2f s\n\n", __func__, nb_treated, Get_T()-startTime);
    }
}
#endif

