#include <stdint.h>
#include <math.h>
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>
#include <assert.h>
#include <pthread.h>
#include "Overlaps.h"
#include "Process_Read.h"
#include "CommandLines.h" 
#include "Overlaps_hamt.h"
#include "ksort.h"

#define HAMT_MAX(x, y) ((x >= y)?(x):(y))  // same
#define HAMT_MIN(x, y) ((x <= y)?(x):(y))  // same
#define HAMT_DIFF(x, y) ((HAMT_MAX((x), (y))) - (HAMT_MIN((x), (y))))  // it's in Correct.h but don't want to include so


KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)
KRADIX_SORT_INIT(ovhamt64, uint64_t, uint64_t, 8)
KRADIX_SORT_INIT(ovhamt32, uint32_t, uint32_t, 4)

#define HAMT_PRIMARY_LABEL 0
#define HAMT_ALTER_LABEL 1
void hamt_ug_util_BFS_markSubgraph_trailing(ma_ug_t *ug_old, ma_ug_t *ug_new, int base_label);


const unsigned char seq_nt4_table[256] = {
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
float cosine_similarity(float *a, float *b, int l){
    float d=0, A=0, B=0;
    for (int i=0; i<l; i++){
        d+= a[i]*b[i];
        A+= a[i]*a[i];
        B+= b[i]*b[i];
    }
    return d/(sqrt(A)*sqrt(B));

}

//////////////////////////////////////////////////////////////////////
//                        debug functions                           //
//////////////////////////////////////////////////////////////////////
// for debug fprintf from kthreads
// (note: printf operates in multithreaded; might crash if race)
void hamt_dbgmsg_init(dbgmsg_t *h){
    h->n = 0;
    h->m = 128;
    h->a = (char*)malloc(h->m * 1);
    assert(h->a);
}
void hamt_dbgmsg_destroy(dbgmsg_t *h){
    free(h->a);
}
void hamt_dbgmsg_reset(dbgmsg_t *h){
    h->n = 0;
}
void hamt_dbgmsg_append(dbgmsg_t *h, char *s, int l){
    while (h->n+l>=h->m){
        h->m = h->m + (h->m>>1);
        h->a = (char*)realloc(h->a, h->m);
        assert(h->a);
    }
    sprintf(h->a+h->n, "%s", s);
    h->n += l;
}
int hamt_dbgmsg_is_empty(dbgmsg_t *h){
    if (h->n==0){return 1;}
    return 0;
}

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
    //    Does not include subgraph info.
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
    q->prev_m = q->m;
}

void queue32_destroy(queue32_t *q){
    free(q->a);
}

void queue32_enqueue(queue32_t *q, uint32_t d){
    if ((q->n+1)>=q->m){
        uint32_t shift = (q->m>>1);
        q->a = (uint32_t*)realloc(q->a, sizeof(uint32_t)*(q->m + shift));
        assert(q->a);
        q->prev_m = q->m;
        q->m = q->m + shift;
        if (q->i_head>q->i_tail){
            // need to shift the content between head and prev_m 
            for (int idx=q->prev_m-1; idx>=q->i_head; idx--){
                q->a[idx+shift] = q->a[idx];
            }
            q->i_head += shift;

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
    if (q->i_head==(q->m-1)){
        q->i_head = 0;
    }else{
        q->i_head++;
    }
    return 1;
}

void queue32_reset(queue32_t *q){
    q->n = 0;
    q->i_head = 0;
    q->i_tail = 0;
    q->prev_m = q->m;
}

int queue32_get_size(queue32_t *q){
    return q->n;
}

int queue32_isempty(queue32_t *q){
    return !(q->n);
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
        if (i==(q->m-1)){
            i = 0;
        }else{
            i++;
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
int stacku32_is_empty(stacku32_t *stack){
    if (stack->n==0) return 1;
    return 0;
}
int stacku32_get_length(stacku32_t *stack){
    return stack->n;
}
int stacku32_peek_last_item(stacku32_t *stack, uint32_t *buf){
    // RETURN
    //     1 if success
    //     0 if fail (stack is empty)
    if (stack->n>0){
        *buf = stack->a[stack->n-1];
        return 1;
    }
    return 0;
}
void stacku32_invert_buffer(stacku32_t *stack){
    if (stack->n<=1) return;
    uint32_t tmp;
    for (int i=0; i<stack->n/2; i++){  // TODO/BUG: is flooring always guaranteed?
        tmp = stack->a[i];
        stack->a[i] = stack->a[stack->n-1-i];
        stack->a[stack->n-1-i] = tmp;
    }
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

void stacku32_copyover(stacku32_t *source, stacku32_t *dest){
    for (int i=0; i<source->n; i++){
        stacku32_push(dest, source->a[i]);
    }
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
int stack32_is_in_stack_givenrange(stacku32_t *stack, uint32_t d, int start, int end){
    // FUNC
    //     check if a given value is in the buffer's segment [start, end)
    // RETURN
    //     1 if yes
    //     0 if no
    if (stack->n==0){return 0;}
    end = end>stack->n? stack->n : end;
    start = start<0? 0 : start;
    for (int i=start; i<end; i++){
        if (stack->a[i]==d){
            return 1;
        }
    }
    return 0;
}

int uint32_buffer_unordered_equal(uint32_t *buf1, int buf1_l, uint32_t *buf2, int buf2_l){
    if (buf1_l!=buf2_l) return 0;
    int ret = 1;
    uint32_t *b1 = (uint32_t*)malloc(buf1_l * sizeof(uint32_t));
    uint32_t *b2 = (uint32_t*)malloc(buf1_l * sizeof(uint32_t));
    memcpy(b1, buf1, buf1_l);
    memcpy(b2, buf2, buf2_l);
    radix_sort_ovhamt32(b1, b1+buf1_l);
    radix_sort_ovhamt32(b2, b2+buf2_l);
    for (int i=0; i<buf1_l; i++){
        if (b1[i]!=b2[i]){
            ret = 0;
            break;
        }
    }
    return ret;
}

int stacku32_unordered_equal(stacku32_t *stack1, stacku32_t *stack2){
    // FUNC
    //     Compare if stack1 and stack2 has the same content
    return uint32_buffer_unordered_equal(stack1->a, stack1->n, stack2->a, stack2->n);
}

int stacku32_index_value(stacku32_t *stack, uint32_t d){
    // FUNC
    //    if a given value is in the buffer, return the index
    //    return -1 otherwise
    if (stack->n==0){return -1;}
    for (int i=0; i<stack->n; i++){
        if (stack->a[i]==d){
            return i;
        }
    }
    return -1;
}


// quick and dirty struct for topo sort
typedef struct {
    uint64_t *a;
    int n, m;  // n is the count, m is the capacity
}stacku64_t;
void stacku64_init(stacku64_t *stack){
    stack->n = 0;
    stack->m = 16;
    stack->a = (uint64_t*)malloc(16*sizeof(uint64_t));
    assert(stack->a);
}
void stacku64_destroy(stacku64_t *stack){
    free(stack->a);
}
void stacku64_push(stacku64_t *stack, uint64_t value){
    if (stack->n+2>stack->m){
        stack->m = stack->m<16? 16 : (stack->m + (stack->m>>1));
        stack->a = (uint64_t*)realloc(stack->a, sizeof(uint64_t)*stack->m);
        assert(stack->a);
    }
    stack->a[stack->n++] = value;
}
int stacku64_pop(stacku64_t *stack, uint64_t *value){
    if (stack->n==0){
        return 0;
    }
    *value = stack->a[stack->n-1];
    stack->n--;
    return 1;
}
void stacku64_reset(stacku64_t *stack){
    stack->n = 0;
}
void stacku64_invert_buffer(stacku64_t *stack){
    if (stack->n<=1) return;
    uint32_t tmp;
    for (int i=0; i<stack->n/2; i++){  // TODO/BUG: is flooring always guaranteed?
        tmp = stack->a[i];
        stack->a[i] = stack->a[stack->n-1-i];
        stack->a[stack->n-1-i] = tmp;
    }
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
    // fprintf(stderr, "- droparc: utg%.6d - utg%.6d\n", (int)(vu>>1)+1, (int)(wu>>1)+1);
}

static inline void hamt_ug_utg_softdel(asg_t *sg, ma_ug_t *ug, uint32_t vu, int label){
    // label unitig vu and its corresponding reads
    uint32_t v;
    if (ug->u.a[vu>>1].len>1000000 && label>0){
        fprintf(stderr, "[W::%s] (tried to mark a very long utg (%.6d l %d) as alt; blocked)\n", __func__, (int)(vu>>1)+1, (int)ug->u.a[vu>>1].len);
        return;
    }
    // fprintf(stderr, "- softdel: utg%.6d (length %d), to label %d\n", (int)(vu>>1)+1, (int)ug->u.a[vu>>1].len, label);
    for (int i=0; i<ug->u.a[vu>>1].n; i++){
        v = (uint32_t) (ug->u.a[vu>>1].a[i]>>33);
        sg->seq[v].c = label;
    }
    ug->u.a[vu>>1].c = label;
    ug->g->seq_vis[vu>>1] = label;  // note seq_vis is actually n_seq*2 long
    
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
    int verbose = 0;

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
            if ((av[i].v>>1)!=(end>>1) && (av[i].v>>1)!=(start>>1)){
                hamt_ug_utg_softdel(sg, ug, av[i].v, new_label);
                if (verbose){
                    fprintf(stderr, "[debug::%s] removed %.6d dir %d\n", __func__, (int)(av[i].v>>1)+1, (int)(av[i].v&1));
                }
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
    // double startTime = Get_T();

    if (ug->utg_coverage){
        free(ug->utg_coverage);
    }
    ug->utg_coverage = (int*)calloc(ug->u.n, sizeof(int));
    uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));
    for (uint32_t i=0; i<ug->u.n; i++){
        ug->utg_coverage[i] = get_ug_coverage(&ug->u.a[i], sg, coverage_cut, sources, ruIndex, primary_flag);
    }
    free(primary_flag);
    // if (VERBOSE>1){
    //     fprintf(stderr, "[M::%s] collected ug coverages.\n", __func__);
    //     fprintf(stderr, "[T::%s] took %0.2f s\n\n", __func__, Get_T()-startTime);
    // }
}
void hamt_destroy_utg_coverage(ma_ug_t *ug){
    if (ug->utg_coverage){
        free(ug->utg_coverage);
        ug->utg_coverage = 0;
    }
}

void hamt_ug_cleanup_arc_by_labels(asg_t *sg, ma_ug_t *ug){
    // (legacy)
    // NOTE
    //     ma_ug_gen_primary only looks at the label when seeding unitig,
    //     so for soft cut we still need to drop some arcs.
    //     Normaly routines should handle arc removal themselves,
    //       this function meant for waterproof/debug.
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
    // if (VERBOSE>1){
    //     fprintf(stderr, "[debug::%s] cleaned up %d arcs\n", __func__, cnt);
    // }
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
    hamt_ug_util_BFS_markSubgraph(ug, flag);  // get subgraph IDs
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

    ma_ug_t *ug_new = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, flag);
    hamt_ug_util_BFS_markSubgraph_trailing(*ug, ug_new, flag);

    hamt_ug_destroy(*ug);
    *ug = ug_new;
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

void asg_get_subgraph_DFSfinishTimes(asg_t *sg, stacku32_t *vertices, stacku32_t *buffer){
    // FUNC
    //     Given a subgraph indicated by vertices (a list of vertex IDs),
    //      perform toposort and store the sorted list of vertex IDs in buffer.
    // NOTE
    //     Won't check base_label (sg can be sg or auxsg).
    int verbose = 0;
    if (stacku32_is_empty(vertices)){
        fprintf(stderr, "[W::%s] tried to toposort an empty subgraph\n", __func__);
        return;
    }
    // stacku32_reset(buffer);
    stacku64_t finishing_time;
    stacku64_init(&finishing_time);
    
    // DFS buffer
    stacku32_t stack;
    uint8_t *color = (uint8_t*)calloc(sg->n_seq*2, 1);

    // init
    stacku32_init(&stack);
    uint32_t vu=0;
    stacku32_peek_last_item(vertices, &vu);
    stacku32_push(&stack, vu);
    color[vu] = 1;

    // DFS
    uint32_t time = 0;
    int is_exhausted;
    uint32_t nv;
    asg_arc_t *av;
    while (stacku32_pop(&stack, &vu)){
        if (verbose) {fprintf(stderr, "[debug::%s] pop %.6d\n", __func__, (int)(vu>>1)+1);}
        nv = asg_arc_n(sg, vu);
        av = asg_arc_a(sg, vu);
        is_exhausted = 1;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (color[av[i].v]==0 && stack32_is_in_stack(vertices, av[i].v)){
                is_exhausted = 0;
                color[av[i].v] = 1;
                stacku32_push(&stack, vu);  // put back vu
                if (verbose) {fprintf(stderr, "[debug::%s]    put back\n", __func__);}
                if (verbose) {fprintf(stderr, "[debug::%s]    push %.6d\n", __func__, (int)(av[i].v>>1)+1);}
                stacku32_push(&stack, av[i].v);
                break;
            }
        }

        if (is_exhausted){
            color[vu] = 2;
            stacku64_push(&finishing_time, ((uint64_t)time)<<32 | vu );
            if (verbose) {fprintf(stderr, "[debug::%s] finished %.6d , t is %d\n", __func__, (int)(vu>>1)+1, (int)(time));}
        }

        if (time<UINT32_MAX) {time++;}
        else{fprintf(stderr, "[W::%s] too many steps\n", __func__);}

        if (stacku32_is_empty(&stack)){  // try to get a new seed
            for (int i=0; i<vertices->n; i++){
                if (color[vertices->a[i]]==0){
                    vu = vertices->a[i];
                    stacku32_push(&stack, vu);
                    color[vu] = 1;
                    if (time<UINT32_MAX) {time++;}
                    else{fprintf(stderr, "[W::%s] too many steps\n", __func__);}
                    if (verbose){
                        fprintf(stderr, "[debug::%s] seed %.6d\n", __func__, (int)(vu>>1)+1);
                    }
                    break;
                }
            }
        }
    }

    // collected relevant finishing times
    uint32_t t;
    int sancheck_counter = 0;
    for (int i=0; i<finishing_time.n; i++){
        for (int j=0; j<vertices->n; j++){
            vu = (uint32_t)finishing_time.a[i];
            t = (uint32_t) (finishing_time.a[i]>>32);
            if (vu==vertices->a[j]){
                buffer->a[j] = t;
                sancheck_counter+=1;
                break;
            }
        }
    }
    assert(sancheck_counter==vertices->n); // all requested vertices shall have a finishing time

    // cleanup
    free(color);
    stacku64_destroy(&finishing_time);
}

void hamt_util_shortest_path(){
    ;
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
    int verbose = 0;

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
    if (verbose>1){
        fprintf(stderr, "[debug::%s] check utg%.6d (shorter) - utg%.6d\n", __func__, (int)(vu_short>>1)+1, (int)(vu_long>>1)+1);
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
                    if (verbose>1){
                        fprintf(stderr, "[debug::%s]     read pair: %.*s - %.*s\n", __func__, 
                                        (int)Get_NAME_LENGTH(R_INF, qn), Get_NAME(R_INF, qn),
                                        (int)Get_NAME_LENGTH(R_INF, qn2), Get_NAME(R_INF, qn2));
                    }
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
    int ret = 0;
    if (cnt_hit==0){ret = -1;}
    if ((float)cnt_hit/cnt_total > ratio ){ret = 1;}

    if (verbose){fprintf(stderr, "[debug::%s] ret is %d\n", __func__, ret);}
    return ret;
}

void hamt_check_diploid_report(ma_ug_t *ug, uint32_t vu1, uint32_t vu2, ma_hit_t_alloc *reverse_sources,
                              float *ratio1, float *ratio2)
{
    // FUNC
    //      debug version of hamt_check_diploid, report stuff instead of return 0 or 1
    // RETURN
    //     -1 if no hit at all
    //     1 if yes
    //     0 if no
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu_short, vu_long;
    uint32_t qn, qn2, tn;
    int nb_targets;
    int found, cnt_hit=0;
    vecu32_t v;
    vecu32_init(&v);
    
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
                    if (!vecu32_is_in_vec(&v, qn2)){
                        vecu32_push(&v, qn2);
                    }
                    break;
                }
            }
        }

        // update stats
        if (found){
            cnt_hit++;
        }
    }
    
    *ratio1 = (float)cnt_hit / ug->u.a[vu_short>>1].n;
    *ratio2 = (float)v.n / ug->u.a[vu_long>>1].n;
    vecu32_destroy(&v);
}

int hamt_check_suspicious_diploid(asg_t *sg, ma_ug_t *ug, uint32_t vu1, uint32_t vu2, float ratio){
    // FUNC
    //     (use case: vu1 and vu2 aren't normal haplotigs - determined by the caller, not checked here.) 
    //     Check OVEC_INFO to see if vu1 and vu2 had some not-that-good overlaps.
    // RETURN
    //     1 if yes
    //     0 if no
    //     (TODO)-1 if overlap state is very asymmetric or other unideal situation
    // TODO
    //     check read direction
    int verbose = 0;

    vecu32_t tig1_reads, tig2_reads, targets;  // note: has direction
    vecu32_init(&tig1_reads);
    vecu32_init(&tig2_reads);
    vecu32_init(&targets);
    uint32_t qn, tn;
    int hits = 0, total;
    ovecinfo_t *h;
    // collected reads (bc need linear search)
    for (int i=0; i<ug->u.a[vu1>>1].n; i++){
        vecu32_push(&tig1_reads, (uint32_t)(ug->u.a[vu1>>1].a[i]>>33));
    }
    for (int i=0; i<ug->u.a[vu2>>1].n; i++){
        vecu32_push(&tig2_reads, (uint32_t)(ug->u.a[vu2>>1].a[i]>>33));
    }
    total = tig1_reads.n<tig2_reads.n? tig1_reads.n : tig2_reads.n;

    int passed=0;
    for (int i=0; i<tig1_reads.n; i++){
        passed = 0;
        qn = tig1_reads.a[i];
        for (int j=0; j<R_INF.OVEC_INF.a[qn].n; j++){
            tn = R_INF.OVEC_INF.a[qn].tn[j];
            if (vecu32_is_in_vec(&tig2_reads, tn)){
                passed = 1;
                if (!vecu32_is_in_vec(&targets, tn)){
                    vecu32_push(&targets, tn);
                }
            }
        }
        if (passed){hits++;}
    }
    int ret;
    if ((float)passed/total>ratio || (float)targets.n/tig2_reads.n>ratio){
        ret = 1;
    }else{
        ret = 0;
    }

    if (verbose){
        fprintf(stderr, "[debug::%s] ret %d, rate %.2f\n", __func__, ret, (float)passed/total);
        fprintf(stderr, "[debug::%s] trace targeted reads:\n", __func__);
        for (int i=0; i<targets.n; i++){
            fprintf(stderr, "[debug::%s]     %.*s\n", __func__, (int)Get_NAME_LENGTH(R_INF, targets.a[i]), Get_NAME(R_INF, targets.a[i]));
        }
    }

    vecu32_destroy(&tig1_reads);
    vecu32_destroy(&tig2_reads);
    vecu32_destroy(&targets);
    return ret;
}
void hamt_check_suspicious_diploid_report(asg_t *sg, ma_ug_t *ug, uint32_t vu1, uint32_t vu2, 
                                          float *ratio1, float *ratio2){
    // FUNC
    //     debug version, report the ratios
    int verbose = 0;

    vecu32_t tig1_reads, tig2_reads, targets;  // note: has direction
    vecu32_init(&tig1_reads);
    vecu32_init(&tig2_reads);
    vecu32_init(&targets);
    uint32_t qn, tn;
    int hits = 0, total;
    ovecinfo_t *h;
    // collected reads (bc need linear search)
    for (int i=0; i<ug->u.a[vu1>>1].n; i++){
        vecu32_push(&tig1_reads, (uint32_t)(ug->u.a[vu1>>1].a[i]>>33));
    }
    for (int i=0; i<ug->u.a[vu2>>1].n; i++){
        vecu32_push(&tig2_reads, (uint32_t)(ug->u.a[vu2>>1].a[i]>>33));
    }
    total = tig1_reads.n<tig2_reads.n? tig1_reads.n : tig2_reads.n;

    int passed=0;
    for (int i=0; i<tig1_reads.n; i++){
        passed = 0;
        qn = tig1_reads.a[i];
        for (int j=0; j<R_INF.OVEC_INF.a[qn].n; j++){
            tn = R_INF.OVEC_INF.a[qn].tn[j];
            if (vecu32_is_in_vec(&tig2_reads, tn)){
                passed = 1;
                if (!vecu32_is_in_vec(&targets, tn)){
                    vecu32_push(&targets, tn);
                }
            }
        }
        if (passed){hits++;}
    }

    *ratio1 = (float)hits / tig1_reads.n;
    *ratio2 = (float)targets.n / tig2_reads.n;

    vecu32_destroy(&tig1_reads);
    vecu32_destroy(&tig2_reads);
    vecu32_destroy(&targets);
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

int hamt_check_if_share_target(asg_t *g, uint32_t v, uint32_t u, int ignore_direction, int base_label){
    uint32_t nv, nu;
    asg_arc_t *av, *au;
    nv = asg_arc_n(g, v);
    av = asg_arc_a(g, v);
    nu = asg_arc_n(g, u);
    au = asg_arc_a(g, u);
    for (int i=0; i<nv; i++){
        if (av[i].del) {continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label) {continue;}
        for (int j=0; j<nu; j++){
            if (au[j].del) {continue;}
            if (base_label>=0 && g->seq_vis[au[j].v>>1]!=base_label) {continue;}
            if (av[i].v==au[j].v){
                return 1;
            }
            if (ignore_direction && ( (av[i].v>>1) == (au[j].v>>1) )){
                return 1;
            }
        }
    }   
    return 0;
}
int hamt_check_if_is_immediate_decedent_of(asg_t *g, uint32_t parent, uint32_t child){
    uint32_t nv = asg_arc_n(g, parent);
    asg_arc_t *av = asg_arc_a(g, parent);
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del) continue;
        if (av[i].v==child) return 1;
    }
    return 0;
}
int hamt_check_if_jump_reachable_simple(asg_t *g, uint32_t v, uint32_t u, uint32_t w, int ban_w,
                                   uint8_t ignore_direction){
    // FUNC
    //     Check if node v can reach node u easily, which is defined as 
    //       the pass only containing one additional node that is NOT w.
    //     Note that dirrect v->u does not count (or we don't need this function).
    // PARAMETER
    //     Note that v and u should always contain the direction bit.
    // TODO?
    //     base_label?
    // RETURN
    //     Number of possible passes (0 if none);
    int ret = 0;
    int verbose = 0;
    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    uint32_t nw;
    asg_arc_t *aw;
    uint8_t ig = !!ignore_direction;

    if (verbose) {fprintf(stderr, "[debug::%s] > check %.6d\n", __func__, (int)(v>>1)+1);}
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del) continue;
        if (ban_w && av[i].v==w) continue;
        if (verbose) {fprintf(stderr, "[debug::%s]  (inter: %.6d)\n", __func__, (int)((av[i].v)>>1)+1);}
        nw = asg_arc_n(g, av[i].v);
        aw = asg_arc_a(g, av[i].v);
        for (uint32_t j=0; j<nw; j++){
            if (aw[j].del) continue;
            if (verbose) {fprintf(stderr, "[debug::%s]  (target: %.6d)\n", __func__, (int)((aw[j].v)>>1)+1);}
            if ((aw[j].v>>ig)==(u>>ig)){
                if (verbose){fprintf(stderr, "[debug::%s]   via %.6d\n",__func__, (int)(aw[j].v>>1)+1);}
                ret++;
            }
        }
    }
    return ret;
}
int hamt_check_if_reachable_within_x_step(ma_ug_t *ug, 
                                          uint32_t start, uint32_t *ends, int number_of_end,
                                          uint8_t *forbidden_buffer, 
                                          int max_step){
    // FUNC
    //     A simple way to check local topology when we don't really need to measure _the _ shortest path.
    //     Check if a path v->...->u exists. It should also: 1) doesn't involve vertex `w` (if ban_w),
    //      and 2) is less than `max_step` steps.
    // PARAMETER
    //     `forbidden_buffer` is an array of length ug->g->n_seq , each entry is either 0 or 1.
    //      0 means the reaching path can contain this unitig, 1 otherwise. 
    // NOTE
    //     This function is for string graph OR unitig graph. If ug is NULL, assume
    //      `sg` is supplied; otherwise will ignore `sg`.
    //     Only ug will be checked for `threshold_max_bp` though.
    // RETURN
    //     0 if no
    //     the minimum step if yes
    int verbose = 0;
    int ret = 0;
    int stepped = 0;
    number_of_end = number_of_end<=0? 1:number_of_end;

    asg_t *sg = ug->g;

    // aux
    stacku32_t stacks[2];
    stacku32_init(&stacks[0]);
    stacku32_init(&stacks[1]);
    int which = 0;  // which stack we are popping stuff from; the other one, i.e. stacks[!which], will receive pushes.
    uint32_t vu=start, nv;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(sg->n_seq*2, 1);  // we care about direction

    // init
    stacku32_push(&stacks[which], vu);
    color[vu] = 1;

    // main loop
    // (switch between two stacks, each step is one 'step')
    while (1){
        while ( stacku32_pop(&stacks[which], &vu) ){
            if ((vu == start) && stepped>0 ){  // ignore paths that loop back
                continue; 
            }
            if (verbose) fprintf(stderr, "[debug::%s] at %.6d\n", __func__, (int)(vu>>1)+1);
            nv = asg_arc_n(sg, vu);
            av = asg_arc_a(sg, vu);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if (color[av[i].v]==0){
                    if (verbose) fprintf(stderr, "[debug::%s]   new encounter: %.6d\n", __func__, (int)(av[i].v>>1)+1);
                    for (int i_end=0; i_end<number_of_end; i_end++){
                        if (av[i].v==(*(ends+i_end))) {  // if you segfault here, check sanity of `number_of_end`
                            if (verbose) fprintf(stderr, "[debug::%s]   FOUND %.6d\n", __func__, (int)(vu>>1)+1);
                            ret = 1;
                            goto finish;
                        }
                    }
                    if (forbidden_buffer && forbidden_buffer[av[i].v>>1]) continue;  // check this only after termination has been checked
                    stacku32_push(&stacks[!which], av[i].v);
                    color[av[i].v] = 1;
                }
            }
            color[vu] = 2;
        }
        stepped++;
        if (stepped>=max_step) goto finish;

        which = !which;
    }
finish:
    free(color);
    stacku32_destroy(&stacks[0]);
    stacku32_destroy(&stacks[1]);
    return ret;
}

int hamt_asgarc_util_isTip(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label);
int hamt_ug_check_if_arc_strongly_bonded(ma_ug_t *ug, uint32_t handle, uint32_t end, 
                                         uint8_t *forbidden_buffer, int max_step){
    // FUNC
    //     "Strongly bonded" is defined as: given that arc vu->wu exists, 
    //       either wu^1 doesn't have targets other than vu^1,
    //       or all the other targets can reach vu or vu's predecessors 
    //       within `max_step` (passing through wu or wu^1 is fine).
    //     This is an attempt to identify the "branching" of convoluted but dense paths
    //       in the bandage layout, without requiring the graph layout.
    //     This check is asymetric! (handle, end) being strong doesn't imply (end^1, handle^1) is strong.
    // NOTE
    //     Please don't supply a huge max_step, there's a little linear search.
    //     Doesn't check if vu->wu exists. You can check vu and wu in vu->uu->wu
    //      using this func if needed.
    // TODO
    //     Check pairs when called for, not the most efficient way right now. Topo sort?
    // RETURN
    //     0 if no
    //     1 if yes
    int verbose = 0;
    int ret = 1;
    asg_t *auxsg = ug->g;
    
    if (verbose) fprintf(stderr, "[debug::%s] handle %.6d dir %d ; end %.6d dir %d \n", __func__,
                                   (int)(handle>>1)+1, (int)(handle&1),
                                   (int)(end>>1)+1, (int)(end&1));
    if (hamt_asgarc_util_countPre(auxsg, end, 0, 0, -1)<=1) {
        if (verbose) fprintf(stderr, "[debug::%s]   end has no backward braching\n", __func__);
        return 1;
    }
    if (hamt_asgarc_util_isTip(auxsg, end, 0, 0, -1)) {
        if (verbose) fprintf(stderr, "[debug::%s]   end is a tip\n", __func__);
        return 1;
    }

    uint32_t nv, vu, wu;
    asg_arc_t *av;

    // a lazy way to collect predecessors of handle
    stacku32_t stacks[3];  // first two are for pushing/popping, while the 3rd collects all. 
    int which = 0;
    int stepped = 0;
    for (int i=0; i<3; i++){stacku32_init(&stacks[i]);}
    vu = handle^1;
    stacku32_push(&stacks[which], vu);
    while (stepped<=max_step){
        while (stacku32_pop(&stacks[which], &vu)){
            if ((vu == (handle^1)) && stepped>0 ){  // ignore paths that loop back
                continue; 
            }
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if (av[i].v==(end^1)){  // don't need to do the search, we have a loop and it connects the `end` node.
                    continue;
                }
                if (!stack32_is_in_stack(&stacks[2], av[i].v)){
                    stacku32_push(&stacks[!which], av[i].v);
                    stacku32_push(&stacks[2], av[i].v);
                }
            }
        }
        which = !which;
        stepped++;
    }
    if (verbose) fprintf(stderr, "[debug::%s]   collected handle's predecessors, count is %d\nthepredecessors:\n", __func__, stacks[2].n);
    if (verbose){
        for (int i=0; i<stacks[2].n; i++){
            fprintf(stderr, "%.6d\n", (int)(stacks[2].a[i]>>1)+1);
        }
    }

    // try each of the backward branches
    nv = asg_arc_n(auxsg, end^1);
    av = asg_arc_a(auxsg, end^1);
    ret = 0;
    for (uint32_t i=0; i<nv; i++){
        if (av[i].del) continue;
        if (av[i].v==(handle^1)) continue;
        if (av[i].v==handle) continue;
        if (verbose)fprintf(stderr, "[debug::%s]   > check branch %.6d\n", __func__, (int)(av[i].v>>1)+1);
        if (hamt_check_if_reachable_within_x_step(ug, av[i].v, stacks[2].a, stacks[2].n, forbidden_buffer, max_step)){
            if (verbose)fprintf(stderr, "[debug::%s]   > passed via branch %.6d\n", __func__, (int)(av[i].v>>1)+1);
            ret = 1;
            goto finish;
        }
    }
    if (verbose)fprintf(stderr, "[debug::%s]   > FAILED all branches\n", __func__);

    

finish:
    for (int i=0; i<3; i++){stacku32_destroy(&stacks[i]);}
    return ret;
   
}
int hamt_ug_check_if_arc_strongly_bonded_bothdir(ma_ug_t *ug, uint32_t handle, uint32_t end, 
                                         uint8_t *forbidden_buffer, int max_step){
    int verbose = 1;
    int ret1 = hamt_ug_check_if_arc_strongly_bonded(ug, handle, end, forbidden_buffer, max_step);
    int ret2 = hamt_ug_check_if_arc_strongly_bonded(ug, end^1, handle^1, forbidden_buffer, max_step);
    if (verbose) {fprintf(stderr, "[debug::%s] handle %.6d end %.6d , 1st dir is_strong=%d, 2nd=%d\n", 
                                    __func__, (int)(handle>>1)+1, (int)(end>>1)+1, ret1, ret2);}
    return (ret1&&ret2);
}

int hamt_ug_arc_del_selfcircle(asg_t *sg, ma_ug_t *ug, uint32_t vu, int base_label){
    // NOTE: caller is responsible for more topo context sanchecks.
    asg_t *auxsg = ug->g;
    // if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)==1 &&
    //     hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==1){
    //         fprintf(stderr, "[W::%s] tried to cut an unconnected circle\n", __func__);
    //         return 0;
    // }
    uint32_t nv;
    asg_arc_t *av;

    hamt_ug_arc_del(sg, ug, vu, vu, 1);
    return 1;
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

int hamt_asgarc_util_get_the_one_target(asg_t *g, uint32_t v, uint32_t *w, int include_del_seq, int include_del_arc, int base_label);
int hamt_asgarc_util_isIsolatedCircle(asg_t *g, uint32_t v, int base_label){
    // FUNC
    //     check if v is an isolated circle
    // RETURN
    //     0 if no
    //     1 if yes
    if (hamt_asgarc_util_countSuc(g, v, 0, 0, base_label)!=1 || hamt_asgarc_util_countPre(g, v, 0, 0, base_label)!=1) return 0;
    uint32_t w;
    hamt_asgarc_util_get_the_one_target(g, v, &w, 0, 0, base_label);
    if (w==v) return 1;
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
    // v has exactly one target (sucessor)
    // SIDE EFFECT
    //     store the target vertex ID in w
    // RETURN
    //     av idx of the target vertex
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

void hamt_asgarc_util_get_the_two_targets(asg_t *g, uint32_t v, uint32_t *w, uint32_t *u,
                                         int include_del_seq, int include_del_arc, int base_label){
    int nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    uint32_t buf[2];
    int idx=0;
    for (i=0; i<(int)nv; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
        buf[idx] = av[i].v;
        idx++;
        if (idx>2){
            fprintf(stderr, "[E::%s] more than two targets, continue anyway\n", __func__);
        }
    }
    assert(idx=2);  // deadly
    *w = buf[0];
    *u = buf[1];
}

int hamt_asgarc_util_has_self_arc(asg_t *g, uint32_t v, int base_label){
    // FUNC
    //     check if v's arcs are just v->v^1 and v^1->v (NOT v->v which is invalid.)
    //     does not included anything marked as deleted.
    // RETURN
    //     0 if not (note that this includes the case where v has no arc)
    //     1 if yes and only self-arc exists (unconnected circular contig)
    //     2 if self-arc and other arcs
    uint32_t tmp1, tmp2;
    int self_arc = 0;
    int other_arc = 0;

    uint32_t nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    for (int i=0; i<nv; i++){
        if (av[i].del) continue;
        if (g->seq[av[i].v>>1].del) continue;
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label) continue;
        
        if (av[i].v==v) self_arc++;
        else other_arc++;
    }
    // 2nd direction
    nv = asg_arc_n(g, v^1);
    av = asg_arc_a(g, v^1);
    for (int i=0; i<nv; i++){
        if (av[i].del) continue;
        if (g->seq[av[i].v>>1].del) continue;
        if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label) continue;
        
        if (av[i].v==(v^1)) self_arc++;
        else other_arc++;
    }

    if (self_arc==0) return 0;
    else{
        if (other_arc==0) return 1;
        else return 2;
    }
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
    //    Contigs that have more than one target on one side with the other side having no targets
    //      will be considered as tips.
    //    Use hamt_asgarc_util_countNoneDanglingTipSuc if you want to rule them out.
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
    //    Only contig with exact 1 target (the query) will be considered as a tip.
    //    Use hamt_asgarc_util_countNoneTipSuc if you want to include multi-target tips.
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
        // if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
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
        if (verbose){fprintf(stderr, "[debug::%s] v-suc\n", __func__);}
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
        fprintf(stderr, "[W::%s] assertion failed, idx is %d\n", __func__, idx);
        return 0;
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
    if (hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(g, w[0], &u1, 0, 0, base_label)==-1){fprintf(stderr, "[%s] shouldn't happen\n", __func__);return 0;}
    if (hamt_asgarc_util_get_the_one_target_ignoreDanglingTip(g, w[1], &u2, 0, 0, base_label)==-1){fprintf(stderr, "[%s] shouldn't happen\n", __func__);return 0;}
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

    uint32_t s, e, e2, sib=0;  // start, end, sibling
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
    if (i==nv){
        fprintf(stderr, "[W::%s] shouldn't happen - failed to get the sibling edge; continue anyway\n", __func__);
        return 0;
    }
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
        idx = 0;
        for (uint32_t i2=0; i2<asg_arc_n(sg, w); i2++){  // find the non-tip target
            if (aw[i2].del){continue;}
            if (hamt_asgarc_util_isTip(sg, aw[i2].v, 0, 0, base_label)){continue;}
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
            ----w0--
          /         \
    .....v0         |
          \        /
           -----u0-----...
          we want to drop v0->w1 and w1->u1^1
    */
   //     only checkout the most simple case, any branching and it'll be ignored
   //     u0's target can be v0.
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
    if (verbose){fprintf(stderr, "[debug::%s]     at ID %.6d\n", __func__, (int)(v0>>1)+1);}
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
        idx_rev = 0;  // uw[0] is the w0 in the topo comment
    }else if (nuw[0]==2 && nuw[1]==1){
        idx_rev = 1;  // uw[1] is the w0
    }else{
        if (verbose){fprintf(stderr, "[debug::%s]    failed 1-and-2\n", __func__);}
        return 0;
    }
    // check topo: the inversion
    if (hamt_asgarc_util_get_the_one_target(sg, uw[idx_rev], &x, 0, 0, base_label)<0){
        fprintf(stderr, "[W::%s]     failed to get target when supposed to.\n", __func__);fflush(stderr);
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

int hamt_asgarc_util_checkSimplePentaBubble(ma_ug_t *ug, asg_t *g, uint32_t v0, uint32_t *buffer, int base_label){
    // FUNC
    //    pops the following structure:
    /*       v2---v5.
           /   \    \ 
       ---v1    v4  v7---
           \     \ /
            v3---v6
        (or the mid link be v3->v4->v5 instead of v2->v4->v6. it's the same.)
        note that:
          v2 and v3 can NOT have other connections
          v5 and v6 CAN have other connections (either direction) 
    */
   //     by removing the 2 links (4 arcs) in path v2-v4-v6 
   // NOTE
   //    buffer is at least length 3.
   //    doesn't tolerate tips.
   //    also assumes only arc can be marked as del, not sequences (like most other funcs in hamt graph cleaning)
   // NOTE2
   //    if ug is not NULL, will assume g is auxsg and check contig length.
   // RETURN
   //    1 if found
   //    0 if none
    int verbose = 0;
    uint32_t w1[2], w2[2], w3[2], w_tmp;
    uint32_t v2, v3, v4, v5, v6, v7, v7_;
    uint32_t nv;
    asg_arc_t *av;
    int idx;
    int max_contig_length = 100000;
    if (verbose){fprintf(stderr, "[debug::%s] at %.6d\n", __func__, (int)(v0>>1)+1);}

    if (hamt_asgarc_util_countSuc(g, v0, 0, 0, base_label)!=2){
        if (verbose){fprintf(stderr, "    v0 failed\n");}
        return 0;
    }
    // try to get v2 and v3
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
        if (verbose) {fprintf(stderr, "[W::%s] abnormal idx (is %d), check soft failed\n", __func__, idx);}
        return 0;
    }
    // assign v2 and v3
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
        if (verbose) {fprintf(stderr, "[W::%s] failed at v3->v6, continue anyway\n", __func__);}
        return 0;
    }
    // if (hamt_asgarc_util_countPre(g, v6, 0, 0, base_label)!=2 || hamt_asgarc_util_countSuc(g, v6, 0, 0, base_label)!=1){
    //     if (verbose){fprintf(stderr, "    v6 topo failed\n");}
    //     return 0;
    // }

    // check v4
    int passed = 0;
    hamt_asgarc_util_get_the_two_targets(g, v2, &w2[0], &w2[1], 0, 0, base_label);
    // is w2[0] v4?
    if (hamt_asgarc_util_countSuc(g, w2[0], 0, 0, base_label)==1 && hamt_asgarc_util_countPre(g, w2[0], 0, 0, base_label)==1){
        hamt_asgarc_util_get_the_one_target(g, w2[0], &w_tmp, 0, 0, base_label);
        if (w_tmp==v6){
            passed = 1;
            v4 = w2[0];
            v5 = w2[1];
        }
    }
    // is w2[1] v4?
    if (hamt_asgarc_util_countSuc(g, w2[1], 0, 0, base_label)==1 && hamt_asgarc_util_countPre(g, w2[1], 0, 0, base_label)==1){
        hamt_asgarc_util_get_the_one_target(g, w2[1], &w_tmp, 0, 0, base_label);
        if (w_tmp==v6){
            if (passed){  // they can't be both v4, abort
                if (verbose) {fprintf(stderr, "[W::%s] topo is weird, abort\n", __func__);}
                return 0;
            }
            passed = 1;
            v4 = w2[1];
            v5 = w2[0];
        }
    }
    if (!passed){return 0;}

    // check if v5 and v6 shares any target (v7)
    if (!hamt_check_if_share_target(g, v5, v6, 0, base_label)){
        if (verbose){fprintf(stderr, "    sink failed\n");}
        return 0;
    }

    // dont need to check v5 and v6's topo
    
    buffer[0] = v2;
    buffer[1] = v4;
    buffer[2] = v6;

    if (ug){  // restrict contig length
        if (ug->u.a[v0>>1].len>max_contig_length ||
            ug->u.a[v2>>1].len>max_contig_length ||
            ug->u.a[v3>>1].len>max_contig_length ||
            ug->u.a[v4>>1].len>max_contig_length ||
            ug->u.a[v5>>1].len>max_contig_length ||
            ug->u.a[v6>>1].len>max_contig_length){
            if (verbose){fprintf(stderr, "    bubble contig length failed\n");}
            return 0;
        }
    }

    return 1;
}

int hamt_ug_check_localSourceSinkPair(ma_ug_t *ug, uint32_t v0, uint32_t w0, 
                                      int *max_length, int *max_coverage,  // side effect
                                      int base_label,
                                      int max_visited)  // limit search span size
{
    // FUNC
    //     Given unitig v0 and w0, check if they are a pair of local source/sink.
    //     Meant for tangle resolving.
    // NOTE
    //    Directions are: v0->(stuff in between)<-w0
    int verbose = 0;
    // if (((v0>>1)+1) == 1374){
    //     verbose = 2;
    // }
    if (verbose){
        fprintf(stderr, "[debug::%s] utg%.6d vs utg%.6d\n", __func__, (int)(v0>>1)+1, (int)(w0>>1)+1);
    }

    asg_t *auxsg = ug->g;
    if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0, base_label)==0 ||
        hamt_asgarc_util_countSuc(auxsg, w0, 0, 0, base_label)==0){
            return 0;
    }

    uint32_t source=v0, sink=w0;
    uint32_t source_nodir=(v0>>1), sink_nodir=(w0>>1);
    uint32_t vu, wu, nv;
    asg_arc_t *av;
    int nb_visited=0, ret=0;
    
    vecu32_t buf_source, buf_sink;
    vecu32_init(&buf_source);
    vecu32_init(&buf_sink);
    
    queue32_t q;
    queue32_init(&q);
    uint8_t *color;
    color = (uint8_t*)calloc(auxsg->n_seq*2, 1);

    // arbitrary check: if the handle hints circular walk, check handle length
    if (source_nodir==sink_nodir){
        if (ug->u.a[source_nodir].len<500000){
            if  (verbose){fprintf(stderr, "[debug::%s]     hinted circular, handle too short\n", __func__);}
            ret = 0;
            goto finish;
        }
    }
    
    // note: IDs in queue have the direction bit; IDs in stack don't

    queue32_enqueue(&q, source);
    color[source] = 1;
    nb_visited = 0;
    // collect all vertices reacheable (undirected) by source without passing though sink (either direction)
    while (queue32_dequeue(&q, &vu)){
        if  (verbose>1){fprintf(stderr, "[debug::%s]     (source) dequeue utg%.6d\n", __func__, (int)(vu>>1)+1);}
        if (color[vu]==2) {continue;}

        if ((vu>>1)!=source_nodir && (vu>>1)!=sink_nodir){
            if ( (ug->u.a[vu>>1].len > (*max_length)) && 
                 (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)>0) && 
                 (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)>0) ){
                if (verbose){
                    fprintf(stderr, "[debug::%s] log a new max length, prv %d now %d (vu %.6d), vu has %d targets and %d pre\n",
                        __func__, *max_length, ug->u.a[vu>>1].len, (int)(vu>>1)+1, 
                        hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label),
                        hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)
                    );
                }
                *max_length = (int)ug->u.a[vu>>1].len;
            }
            if (ug->utg_coverage[vu>>1] > (*max_coverage)){
                *max_coverage = ug->utg_coverage[vu>>1];
            }
        }
        
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label){continue;}
            wu = av[i].v;
            // forbid branching to outside of the persumed tangle span
            if (wu==source || wu==sink){  // note yes check direction
                if (source_nodir==sink_nodir){
                    ;
                }else{  // not hinted circular, thus forbidden
                    if  (verbose){fprintf(stderr, "[debug::%s]     hit handle (source)\n", __func__);}
                    ret = 0;
                    goto finish;
                }
            }
            if ((wu>>1)==sink_nodir || (wu>>1)==source_nodir){continue;}
            if (color[wu]==0){
                if  (verbose>1){fprintf(stderr, "[debug::%s]     (source)   push utg%.6d\n", __func__, (int)(wu>>1)+1);}
                queue32_enqueue(&q, wu);
                color[wu] = 1;
                if (color[wu^1]==0){
                    queue32_enqueue(&q, wu^1);
                    color[wu^1] = 1;
                }                
                nb_visited++;
            }
        }
        vecu32_push(&buf_source, vu>>1);
        color[vu] = 2;

        if (nb_visited>max_visited){  // search span too big, assume failure
            if  (verbose){fprintf(stderr, "[debug::%s]     search span too big (source)\n", __func__);}
            ret = 0;
            goto finish;
        }
    }

    queue32_reset(&q);
    queue32_enqueue(&q, sink);
    memset(color, 0, auxsg->n_seq*2);
    color[sink] = 1;
    nb_visited = 0;
    // collect all vertices reacheable by sink without passing though sink (either direction)
    while (queue32_dequeue(&q, &vu)){
        if  (verbose>1){fprintf(stderr, "[debug::%s]     (sink) dequeue utg%.6d\n", __func__, (int)(vu>>1)+1);}
        if (color[vu]==2) {continue;}
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label){continue;}
            wu = av[i].v;
            // forbid branching to outside of the persumed tangle span
            if (wu==source || wu==sink){  // note yes check direction
                if (source_nodir==sink_nodir){
                    ;
                }else{  // not hinted circular, thus forbidden
                    if  (verbose){fprintf(stderr, "[debug::%s]     hit handle (sink)\n", __func__);}
                    ret = 0;
                    goto finish;
                }
            }
            if ((wu>>1)==sink_nodir || (wu>>1)==source_nodir){continue;}
            if (color[wu]==0){
                queue32_enqueue(&q, wu);
                color[wu] = 1;
                if (color[wu^1]==0){
                    queue32_enqueue(&q, wu^1);
                    color[wu^1] = 1;
                }
                if  (verbose>1){fprintf(stderr, "[debug::%s]     (sink)   push utg%.6d\n", __func__, (int)(wu>>1)+1);}
                nb_visited++;
            }
        }
        vecu32_push(&buf_sink, vu>>1);
        color[vu] = 2;

        if (nb_visited>max_visited){  // search span too big, assume failure
        if  (verbose){fprintf(stderr, "[debug::%s]     search span too big (sink)\n", __func__);}
            ret = 0;
            goto finish;
        }
    }

    // check if the two collections are identical /*(in opposite direction)*/
    if (buf_source.n!=buf_sink.n){
        if  (verbose){fprintf(stderr, "[debug::%s]     buf unequal size (%d, %d)\n", __func__, buf_source.n, buf_sink.n);}
        ret = 0;
        goto finish;
    }
    for (int i=0; i<buf_source.n; i++){
        if ((buf_source.a[i])==source_nodir || (buf_source.a[i])==sink_nodir){continue;}
        if (!vecu32_is_in_vec(&buf_sink, buf_source.a[i])){
            if  (verbose){fprintf(stderr, "[debug::%s]     source vertex utg%.6d not in sink search\n", __func__, (int)(buf_source.a[i])+1);}
            ret = 0;
            goto finish;
        }
    }
    for (int i=0; i<buf_sink.n; i++){
        if ((buf_source.a[i])==source_nodir || (buf_source.a[i])==sink_nodir){continue;}
        if (!vecu32_is_in_vec(&buf_source, buf_sink.a[i])){
            if  (verbose){fprintf(stderr, "[debug::%s]     sink vertex utg%.6d not in source search\n", __func__, (int)(buf_sink.a[i])+1);}
            goto finish;
        }
    }
    ret = 1;
    if (verbose){fprintf(stderr, "[debug::%s]     success\n", __func__);}

finish:
    vecu32_destroy(&buf_source);
    vecu32_destroy(&buf_sink);
    queue32_destroy(&q);
    free(color);
    return ret;
}


int hamt_ug_util_popSimpleBiBubbleChain(asg_t *sg, ma_ug_t *ug, uint32_t v0, 
                                        int base_label){
    // TODO
    //     This was written before implementing complex bubble popping.
    //     If we don't care about recovering the bubble edge (doing this; 
    //       if recovered it's still chimeric haplotype anyway),
    //       this function can be safely removed. 
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
    int max_contig_length = 100000;

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
        if (av[i].v==v0){return 0;}  // prevent circling back to the handle
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
        if (w_tmp==v0){return nb_cut;}  // prevent circling back to handle
        if (hamt_asgarc_util_countPre(auxsg, w_tmp, 0, 0, base_label)!=2){
            return nb_cut;
        }
        if (hamt_asgarc_util_countSuc(auxsg, w_tmp, 0, 0, base_label)>2){
            return nb_cut;
        }
        if (verbose){fprintf(stderr, "    checkpoint 1\n");}

        if (ug->u.a[w1[0]>>1].len>max_contig_length ||
            ug->u.a[w1[1]>>1].len>max_contig_length ||
            ug->u.a[w_tmp>>1].len>max_contig_length ){
                return nb_cut;
        }

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
            if (av[i].v==v0){return nb_cut;}  // prevent circling back to handle
            if (ug->u.a[av[i].v>>1].len>max_contig_length){
                return nb_cut;
            }

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

        // (restrict length)
        if (ug->u.a[w2[0]>>1].len>max_contig_length ||
            ug->u.a[w2[1]>>1].len>max_contig_length){
                return nb_cut;
        }

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

int hamt_ug_popSpecialCase1(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label){
    /*
    Treat the following topo:
              --v2--v3---v4-
            /        |    \
     -----v0        |      v6-----
           \        |     /
            --v1---v3^1--v5-
    ... where v1 and v4 connect to v3's one side, and v2 & v5 connect to the other. 
    This cannot be identified by the implemented complex bubble check. (Dec 11 2020, r9) 
    */
    // NOTE
    //    path is arbitrarily chosen.
    //    D5196-T2 pctg #141
    asg_t *auxsg = ug->g;
    int nb_treated = 0;
    uint32_t v1, v2, v3, v4, v5, v6;
    uint32_t v_tmp1, v_tmp2;
    for (uint32_t v0=0; v0<auxsg->n_seq*2; v0++){
        if (hamt_asgarc_util_countSuc(auxsg, v0, 0, 0, base_label)!=2){continue;}
        // get v1 and v2
        hamt_asgarc_util_get_the_two_targets(auxsg, v0, &v1, &v2, 0, 0, base_label);
        // v1 and v2 check:
        if (hamt_asgarc_util_countPre(auxsg, v1, 0, 0, base_label)!=1 || 
            hamt_asgarc_util_countPre(auxsg, v2, 0, 0, base_label)!=1){continue;}  // backward branching
        if (hamt_asgarc_util_countSuc(auxsg, v1, 0, 0, base_label)!=1 || 
            hamt_asgarc_util_countSuc(auxsg, v2, 0, 0, base_label)!=1){continue;}  // forward branching
        hamt_asgarc_util_get_the_one_target(auxsg, v1, &v_tmp1, 0, 0, base_label);
        hamt_asgarc_util_get_the_one_target(auxsg, v2, &v_tmp2, 0, 0, base_label);
        if (v_tmp1 != (v_tmp2^1)){continue;}
        v3 = v_tmp1;
        if ((v3>>1)==(v0>>1) || (v3>>1)==(v1>>1) || (v3>>1)==(v2>>1)){continue;}
        // v3 check
        if (hamt_asgarc_util_countPre(auxsg, v3, 0, 0, base_label)!=2){continue;}
        if (hamt_asgarc_util_countSuc(auxsg, v3, 0, 0, base_label)!=2){continue;}
        // get v4 and v5
        hamt_asgarc_util_get_the_two_targets(auxsg, v3, &v_tmp1, &v_tmp2, 0, 0, base_label);
        if ((v_tmp1>>1)==(v1>>1) || (v_tmp1>>1)==(v2>>1)){
            v4 = v_tmp2;
        }else{
            v4 = v_tmp1;
        }
        hamt_asgarc_util_get_the_two_targets(auxsg, v3^1, &v_tmp1, &v_tmp2, 0, 0, base_label);
        if ((v_tmp1>>1)==(v1>>1) || (v_tmp1>>1)==(v2>>1)){
            v5 = v_tmp2;
        }else{
            v5 = v_tmp1;
        }
        // v4 and v5 check
        if ((v4>>1)==(v5>>1)){continue;}
        if (hamt_asgarc_util_countPre(auxsg, v4, 0, 0, base_label)!=1 || 
            hamt_asgarc_util_countPre(auxsg, v5, 0, 0, base_label)!=1){continue;}  // backward branching
        if (hamt_asgarc_util_countSuc(auxsg, v4, 0, 0, base_label)!=1 || 
            hamt_asgarc_util_countSuc(auxsg, v5, 0, 0, base_label)!=1){continue;}  // forward branching
        hamt_asgarc_util_get_the_one_target(auxsg, v1, &v_tmp1, 0, 0, base_label);
        hamt_asgarc_util_get_the_one_target(auxsg, v2, &v_tmp2, 0, 0, base_label);
        if (v_tmp1 != v_tmp2){continue;}
        v6 = v_tmp1;
        if (/*(v6>>1)==(v0>>1) ||*/ (v6>>1)==(v1>>1) || (v6>>1)==(v2>>1) || 
            (v6>>1)==(v3>>1) || (v6>>1)==(v4>>1) || (v6>>1)==(v5>>1)){continue;}
        if (hamt_asgarc_util_countPre(auxsg, v6, 0, 0, base_label)){continue;}

        // chose path
        hamt_ug_arc_del(sg, ug, v0, v1, 1);
        hamt_ug_arc_del(sg, ug, v1, v3, 1);
        hamt_ug_utg_softdel(sg, ug, v1, 1);

        hamt_ug_arc_del(sg, ug, v3, v4, 1);
        hamt_ug_arc_del(sg, ug, v4, v6, 1);
        hamt_ug_utg_softdel(sg, ug, v4, 1);
        nb_treated++;
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d spots\n", __func__, nb_treated);
    }
   return nb_treated;

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
        if (verbose){
            fprintf(stderr, "[debug::%s] seed utg%.6d\n", __func__, (int)(v0>>1)+1);
        }
        if (base_label>=0 && sg->seq_vis[v0>>1]!=base_label){continue;}
        tmp = 0;
        queue32_reset(q);
        queue32_enqueue(q, v0);
        // queue32_enqueue(q, v0^1);
        color[v0] = 1;
        color[v0^1] = 1;
        tmp++;
        subg_labels[v0] = label;
        subg_labels[v0^1] = label;
        while(queue32_dequeue(q, &v)){
            if (color[v]==2){
                if (verbose) {fprintf(stderr, "[debug::%s]    skip utg%.6d bcs color\n", __func__, (int)(v>>1)+1);}
                continue;
            }
            if (verbose) {fprintf(stderr, "[debug::%s]    dequeue utg%.6d, n %d, m %d, head %d, tail %d\n", __func__, (int)(v>>1)+1, (int)q->n, (int)q->m, (int)q->i_head, (int)q->i_tail);}
            // 1st direction
            if (verbose) {fprintf(stderr, "[debug::%s]    check 1st dir\n", __func__);}
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del){continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(q, w);
                    color[w] = 1;
                    color[w^1] = 1;
                    tmp++;
                    subg_labels[w] = label;
                    subg_labels[w^1] = label;
                    if (verbose) {fprintf(stderr, "[debug::%s]    push utg%.6d, n %d, m %d, head %d, tail %d\n", __func__, (int)(w>>1)+1, (int)q->n, (int)q->m, (int)q->i_head, (int)q->i_tail);}
                }
            }
            // 2nd direction
            if (verbose) {fprintf(stderr, "[debug::%s]    check 2nd dir\n", __func__);}
            nv = asg_arc_n(sg, v^1);
            av = asg_arc_a(sg, v^1);
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del){continue;}
                if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                w = av[i].v;
                if (color[w]==0){
                    queue32_enqueue(q, w);
                    color[w] = 1;
                    color[w^1] = 1;
                    tmp++;
                    subg_labels[w] = label;
                    subg_labels[w^1] = label;
                    if (verbose) {fprintf(stderr, "[debug::%s]    push utg%.6d, n %d, m %d, head %d, tail %d\n", __func__, (int)(w>>1)+1, (int)q->n, (int)q->m, (int)q->i_head, (int)q->i_tail);}
                }
            }
            color[v] = 2;
            color[v^1] = 2;
        }
        label++;
        if (verbose){
            fprintf(stderr, "[debug::%s] subgraph #%d, %d reads/tigs.\n", __func__, label, tmp);
        }
    }

    if (!q0){queue32_destroy(&q_alloc);}
    if (!color0){free(color);}
    return label;
}

int hamt_ug_util_BFS_markSubgraph(ma_ug_t *ug, int base_label){
    // FUNC
    //     Traverse the ug and give each unitig a subgraph ID (arbitrary; 
    //       stable in terms of assembly runs, NOT stable between ug generations).
    // RET
    //     The largest subgraph ID.
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
    return max_label;
}

void hamt_ug_util_BFS_markSubgraph_trailing(ma_ug_t *ug_old, ma_ug_t *ug_new, int base_label){
    // FUNC
    //     Init/update a list of subgraphs IDs for each unitig.
    //     Intention: it might be interesting to know if a few seperated circular contigs/subgraphs
    //                were initially linked together by just looking at the contig names.
    //     For the primary graph, do nothing if the connectivity of a subgraph did not change,
    //                            push a new subID otherwise 
    //     For example, if subgraph #4 is splited into two, then #4's unitigs will be updated with
    //       an additional ID of #1 or #2; the final subgraph ID will be #4.1 or #4.2.
    // self NOTE
    //     ug->u.n == ug->g->n_seq == number of unitigs
    int verbose = 0;
    int *old_labels = 0;
    int *new_labels = 0;
    int old_max_label = 0;
    int is_first_time = 0;
    asg_t *auxsg_old = ug_old->g;
    asg_t *auxsg_new = ug_new->g;

    if (!R_INF.subg_label_trail){  // first time
        is_first_time = 1;
        if (verbose) fprintf(stderr, "[debug::%s] init subg_label_trail\n", __func__);
        R_INF.subg_label_trail = (ma_utg_subg_labels_v*)calloc(1, sizeof(ma_utg_subg_labels_v));
        assert(R_INF.subg_label_trail);
        // note: start/end read of a circular unitig might not be stable, so here we have to have a buffer of length nb_reads,
        //       instead of length nb_unitig. 
        //       Not the best way i suppose. // TODO
        R_INF.subg_label_trail->m = R_INF.total_reads; 
        R_INF.subg_label_trail->n = R_INF.total_reads;
        R_INF.subg_label_trail->a = (ma_utg_subglabels_t*)calloc(R_INF.subg_label_trail->m, sizeof(ma_utg_subglabels_t));
        for (int i=0; i<R_INF.total_reads; i++){
            R_INF.subg_label_trail->a[i].m = 16;
            R_INF.subg_label_trail->a[i].a = (int*)malloc(sizeof(int) * 16);
            for (int j=0; j<16; j++){
                R_INF.subg_label_trail->a[i].a[j] = -1;
            }
        }
    }

    // store old subg labels (to be referred by 1st read's ID)
    if (verbose) fprintf(stderr, "[debug::%s] backup previous subg labels\n", __func__);
    old_labels = (int*)calloc(R_INF.total_reads, sizeof(int));  // use seqID as stable unitig ID; note: this works because we only merge unitigs and never break them up
    new_labels = (int*)calloc(R_INF.total_reads, sizeof(int)); 
    for (int i=0; i<R_INF.total_reads; i++){old_labels[i] = -1;}
    for (int i=0; i<R_INF.total_reads; i++){new_labels[i] = -1;}

    for (uint32_t i=0; i<ug_old->u.n; i++){  // note: number of unitig, not 2*nb_unitig
        if (verbose>1) {
            fprintf(stderr, "[debug::%s]     utg%.6d, subg label was %d\n", __func__, (int)(i)+1, (int)ug_old->u.a[i].subg_label);
            fprintf(stderr, "[debug::%s]        (read ID: %d and %d)\n", __func__, (int)(ug_old->u.a[i].a[0]>>33),
                                                                (int)(ug_old->u.a[i].a[ug_old->u.a[i].n-1]>>33));
            // fprintf(stderr, "[debug::%s]         - %.*s\n", __func__, (int)Get_NAME_LENGTH((R_INF), (ug_old->u.a[i].a[0]>>33)),
            //                                                 Get_NAME((R_INF), (ug_old->u.a[i].a[0]>>33)));
            // fprintf(stderr, "[debug::%s]         - %.*s\n", __func__, (int)Get_NAME_LENGTH((R_INF), (ug_old->u.a[i].a[ug_old->u.a[i].n-1]>>33)),
            //                                                 Get_NAME((R_INF), (ug_old->u.a[i].a[ug_old->u.a[i].n-1]>>33)));
            
        }
        // old_labels[ug_old->u.a[i].a[0]>>33] = ug_old->u.a[i].subg_label;
        // old_labels[ug_old->u.a[i].a[ug_old->u.a[i].n-1]>>33] = ug_old->u.a[i].subg_label;
        for (int j=0; j<ug_old->u.a[i].n; j++){
            old_labels[ug_old->u.a[i].a[j]>>33] = ug_old->u.a[i].subg_label;
        }
        
        if (ug_old->u.a[i].subg_label>old_max_label) old_max_label = ug_old->u.a[i].subg_label;
    }
    for (uint32_t i=0; i<ug_new->u.n; i++){
        if (verbose>1) {
            fprintf(stderr, "[debug::%s]     utg%.6d, subg label is %d\n", __func__, (int)(i)+1, (int)ug_new->u.a[i].subg_label);
            fprintf(stderr, "[debug::%s]        (read ID: %d and %d)\n", __func__, (int)(ug_new->u.a[i].a[0]>>33),
                                                                (int)(ug_new->u.a[i].a[ug_new->u.a[i].n-1]>>33));
            // fprintf(stderr, "[debug::%s]         - %.*s\n", __func__, (int)Get_NAME_LENGTH((R_INF), (ug_new->u.a[i].a[0]>>33)),
            //                                                 Get_NAME((R_INF), (ug_new->u.a[i].a[0]>>33)));
            // fprintf(stderr, "[debug::%s]         - %.*s\n", __func__, (int)Get_NAME_LENGTH((R_INF), (ug_new->u.a[i].a[ug_new->u.a[i].n-1]>>33)),
            //                                                 Get_NAME((R_INF), (ug_new->u.a[i].a[ug_new->u.a[i].n-1]>>33)));
        }
        // new_labels[ug_new->u.a[i].a[0]>>33] = ug_new->u.a[i].subg_label;
        // new_labels[ug_new->u.a[i].a[ug_new->u.a[i].n-1]>>33] = ug_new->u.a[i].subg_label;
        for (int j=0; j<ug_new->u.a[i].n; j++){
            new_labels[ug_new->u.a[i].a[j]>>33] = ug_new->u.a[i].subg_label;
        }
    }

    // (note: assuming ug_new's subg labels are already collected)
    int unchanged, base_subgID, subgID;
    stacku32_t *IDs = (stacku32_t*)malloc(sizeof(stacku32_t));
    stacku32_init(IDs);

    // init trailing IDs if 1st time
    if (is_first_time){
        if (verbose) fprintf(stderr, "[debug::%s] init trailing subg IDs\n", __func__);
        for (int i_utg=0; i_utg<ug_old->u.n; i_utg++){
            subgID = ug_old->u.a[i_utg].subg_label;
            if (verbose>1) fprintf(stderr, "[debug::%s]     unitig %.6d\n", __func__, (int)(i_utg)+1);

            ma_utg_subglabels_t *p;
            for (int i_read=0; i_read<ug_old->u.a[i_utg].n; i_read++){
                p = &R_INF.subg_label_trail->a[ug_old->u.a[i_utg].a[i_read]>>33];
                p->a[0] = subgID;
                p->n++;
            }

            if (verbose) fprintf(stderr, "[debug::%s]       (ok, subg ID %d)\n", __func__, (int)subgID);
        }
    }

    // check each old subgraph and see if every unitig still share a same ID
    {
        if (verbose) fprintf(stderr, "[debug::%s] check each (previous) subgraphs...\n", __func__);
        for (int i_subg=0; i_subg<old_max_label+1; i_subg++){
            if (verbose>1) fprintf(stderr, "[debug::%s]   > subg %d\n", __func__, i_subg);
            unchanged = 1;
            base_subgID = -1;
            stacku32_reset(IDs);
            for (int i_utg=0; i_utg<ug_old->u.n; i_utg++){
                // if (verbose>1) fprintf(stderr, "[debug::%s]    at utg%.6d , read ID %d\n", __func__, (int)i_utg+1, (int)(ug_old->u.a[i_utg].a[0]>>33));
                if (ug_old->u.a[i_utg].subg_label!=i_subg) continue;

                // checking a giving subgraph (old)
                for (int i_read=0; i_read<ug_old->u.a[i_utg].n; i_read++){
                    subgID = new_labels[ug_old->u.a[i_utg].a[i_read]>>33];  // new subg ID
                    // if (verbose>1) fprintf(stderr, "[debug::%s]     utg %.6d, a new subg ID is %d\n", __func__, i_utg+1, subgID);
                    if (base_subgID==-1) base_subgID = subgID;
                    else{
                        if (base_subgID!=subgID){
                            unchanged = 0;
                            if (!stack32_is_in_stack(IDs, (uint32_t)subgID)){
                                stacku32_push(IDs, (uint32_t)subgID);
                            }
                        }
                    }
                }
            }

            if (!unchanged){  // an old subgraph has been splited, push new trailing IDs for the reads
                if (verbose) {fprintf(stderr, "[debug::%s] subgraph %d has been splitted; going to udpate trailing IDs\n", __func__, i_subg);}
                stacku32_push(IDs, (uint32_t)base_subgID);
                int newID;
                for (int i_utg=0; i_utg<ug_old->u.n; i_utg++){
                    if (i_subg!=ug_old->u.a[i_utg].subg_label) continue;

                    // push to the trail of subgIDs for all reads
                    for (int i_read=0; i_read<ug_old->u.a[i_utg].n; i_read++){
                        subgID = new_labels[ug_old->u.a[i_utg].a[i_read]>>33];  // new subg ID
                        newID = stacku32_index_value(IDs, subgID);  // the index after subgraph split
                        ma_utg_subglabels_t *p = &R_INF.subg_label_trail->a[ug_old->u.a[i_utg].a[i_read]>>33];
                        if (p->n==p->m){  // (expand buffer)
                            p->m = p->m + (p->m>>1);
                            p->a = (int*)realloc(p->a, sizeof(int)*p->m);
                        }
                        p->a[p->n] = newID;
                        p->n++;
                    }
                }
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s] subgraph %d unchanged\n", __func__, i_subg);}
            }
        }
    }
    
    free(old_labels);
    free(new_labels);
    stacku32_destroy(IDs);
    free(IDs);
}

int hamt_debug_ug_random_cut_arcs(asg_t *sg, ma_ug_t *ug, int nb_cut){
    // FUNC
    //     Arbitrarily cut nb_cut arcs. For debugging the subgraph trailing ID implementation.
    // RET
    //     1 if ok
    //     0 if no more arcs remain
    asg_t *auxsg = ug->g;
    int cnt = 0;
    uint32_t vu, wu;
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del) {continue;}
        vu = (uint32_t) (auxsg->arc[i].ul>>32);
        wu = auxsg->arc[i].v;
        if (auxsg->seq_vis[vu>>1]!=auxsg->seq_vis[wu>>1]){
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            cnt++;
            if (cnt==nb_cut) break;
        }
    }
    if (cnt<nb_cut) return 0;
    else return 1;

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

int hamt_ug_pop_unevenInvertBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label){
    // Like uneven circle, the invert edge of the invert bubble is short and has coverage diff
    //   with the other edges. Check these two aspects and drop if both checks.
    // NOTE
    //     Assumes coverage has been collected.
    // SPECIAL NOTE
    //     This was not a prominent observation (at least in the sheep) before the 
    //        asg-way-too-complex containment read fix; after that, a bunch of them appeared.
    asg_t *auxsg = ug->g;
    int verbose = 0;

    uint32_t vu, wu, uu;
    float vuc, wuc, uuc;  // coverage
    int nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_checkSimpleInvertBubble(auxsg, vu, &wu, &uu, base_label)){
            if (ug->u.a[wu>>1].len>100000) continue;
            if ( (ug->u.a[uu>>1].len<ug->u.a[wu>>1].len) ||
                  ug->u.a[vu>>1].len<ug->u.a[wu>>1].len) continue;  // wu shall be the shortest one

            vuc = ug->utg_coverage[vu>>1];
            wuc = ug->utg_coverage[wu>>1];
            uuc = ug->utg_coverage[uu>>1];
            if( (wuc/vuc)>0.5 || (wuc/uuc)>0.5) continue; 
            
            // cut
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            hamt_ug_arc_del(sg, ug, wu, uu^1, 1);
            if (base_label>=0) hamt_ug_utg_softdel(sg, ug, wu, alt_label);
            nb_cut++;
            if (verbose){
                fprintf(stderr, "[debug::%s] pop: utg%.6d - utg%.6d - utg%.6d\n", __func__,
                                    (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
            }
        }
    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
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
        }
    }
    free(a.a);

}






///////////////////////////////////////////////////////////////////////////////////////
//                  higher-level routines (unitig graph)                             //
///////////////////////////////////////////////////////////////////////////////////////

int hamt_asgarc_count_leads_to_howmany(asg_t *sg, uint32_t vu0, int threshold_max){
    // FUNC
    //     Given an edge v, count how many edges are _connected_.
    //     Note that this function doesn't care about the bidirection.
    //     If threshold_max is given(i.e. >0), will terminate when counter reaches it.
    // RETURN 
    //     number of vertices reachable.
    int ret = 0;
    queue32_t q;
    queue32_init(&q);
    uint32_t vu=vu0, nv;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(sg->n_seq, 1);
    int is_initial = 1;

    queue32_enqueue(&q, vu0);
    color[vu>>1] = 1;

    while (queue32_dequeue(&q, &vu)){
        if (((vu>>1) == (vu0>>1)) ){  // 1st round special case
            if (is_initial) {is_initial = 0;}
            else continue;
        }

        nv = asg_arc_n(sg, vu);
        av = asg_arc_a(sg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            if (color[av[i].v>>1]==0){
                queue32_enqueue(&q, av[i].v);
                queue32_enqueue(&q, av[i].v^1);
                color[av[i].v>>1] = 1;
                ret++;
                if (threshold_max>0 && ret>=threshold_max){
                    goto finish;
                }
            }
        }

        color[vu>>1] = 2;
    }
finish:
    free(color);
    queue32_destroy(&q);
    return ret;
}


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

int hamt_ug_covcut_falseVertexLoop(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    /*         vu
             /  |
    -------wu   |
            \   |
             vu^1
    given that vu's coverage is not much higher than wu
    */
    int verbose = 0;
    asg_t *auxsg = ug->g;
    int nb_cut = 0;
    
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if ( hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=1 || 
             hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=1 ) {
            if (verbose){fprintf(stderr, "[debug::%s] tried %.6d , failed topo\n", __func__, (int)(vu>>1)+1);}
                continue;
        }

        uint32_t wu1, wu2;
        hamt_asgarc_util_get_the_one_target(auxsg, vu, &wu1, 0, 0, base_label);
        hamt_asgarc_util_get_the_one_target(auxsg, vu^1, &wu2, 0, 0, base_label);
        if (wu1!=wu2) {
            if (verbose){
                fprintf(stderr, "[debug::%s] tried %.6d , failed target (wu1 %.6d wu2 %.6d )\n", 
                        __func__, (int)(vu>>1)+1, (int)(wu1>>1)+1, (int)(wu2>>1)+1);
            }
            continue; 
        }

        if (ug->utg_coverage[vu>>1]<(ug->utg_coverage[wu1>>1]*1.5)){
            hamt_ug_arc_del(sg, ug, vu, wu1, 1);
            hamt_ug_arc_del(sg, ug, vu^1, wu1, 1);
            if (verbose){fprintf(stderr, "[debug::%s] cut for %.6d \n", __func__, (int)(vu>>1)+1);}
            nb_cut++;
        }else{
            if (verbose){fprintf(stderr, "[debug::%s] tried %.6d , failed coverage\n", __func__, (int)(vu>>1)+1);}
        }
    }
    return nb_cut;

}

void hamt_asgarc_ugCovCutDFSCircle(asg_t *sg, ma_ug_t *ug, int base_label)
{
    // NOTE
    //     Do not cut if the non-circular part only as one or a few unitigs

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
            // check the non-circular part, if it has only a few edges, do not cut
            if (hamt_asgarc_count_leads_to_howmany(auxsg, v^1, 10)<10){
                ;
            }else{  // cut
                hamt_ug_arc_del(sg, ug, v, w, 1);
                nb_cut++;
                if (verbose){
                    fprintf(stderr, "[debug::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);
                }
            }
        }
    }
    
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d links.\n", __func__, nb_cut);
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
    uint32_t /*covs[20],*/ cov_min, cov_max, cov_v, cov_tmp;

    stacku32_t covs_s;
    stacku32_init(&covs_s);
    uint32_t *covs = covs_s.a;

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
            stacku32_reset(&covs_s);
            for (int i=0; i<nw; i++){
                if (aw[i].del){continue;}
                if (base_label>=0 && auxsg->seq_vis[aw[i].v>>1]!=base_label) {continue;}
                if ((aw[i].v>>1)==(v>>1)){continue;}  // skip the handle
                if (hamt_asgarc_util_isTip(auxsg, aw[i].v, 0, 0, base_label)){continue;}  // don't use tip's coverage
                stacku32_push(&covs_s, (uint32_t)ug->utg_coverage[aw[i].v>>1]);  // covs[idx] = (uint32_t)ug->utg_coverage[aw[i].v>>1]; idx++;                
            }
            stacku32_push(&covs_s, (uint32_t)ug->utg_coverage[w>>1]);  // covs[idx++] = (uint32_t)ug->utg_coverage[w>>1]; idx++;
            idx = covs_s.n;
            covs = covs_s.a;

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

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    stacku32_destroy(&covs_s);
    fprintf(stderr, "[M::%s] cut %d.\n", __func__, nb_cut);
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
    }
}

int hamt_ugasg_cut_shortTips(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    // NOTE: won't treat multi-target tips
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

        #if 0
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
                if (verbose){
                    fprintf(stderr, "[debug::%s] cut tip: utg%.6d \n", __func__, (int)(vu>>1)+1);
                }
            }
        }
        // if (!is_hard_drop && !need_to_spare){
        //     hamt_ug_utg_softdel(sg, ug, vu, alt_label);
        // }
        #endif
        
        // sancheck+cut: don't touch the tip if it has more than one outgoing arc,
        //               regardless of what the targets look like.
        //               We might want to check up on these tips later or something.
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)>1) continue;
        hamt_asgarc_util_get_the_one_target(auxsg, vu, &wu, 0, 0, base_label);
        hamt_ug_arc_del(sg, ug, vu, wu, 1);


    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d unitig tigs.\n", __func__, nb_cut);
    }
    return nb_cut;
}

int hamt_ug_cut_shortTips_arbitrary(asg_t *sg, ma_ug_t *ug, int max_length, int base_label){
    // FUNC
    //     cut simple tips (1 arc) shorter than max_length, regardless of whether their targets have no non-tip linkage
    // NOTE
    //     Will drop multi-target tip. Ideally don't use this before calling hamt_ug_rescue_bifurTip.
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
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<1){continue;}
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

int hamt_ug_check_bifurTip(ma_ug_t *ug, uint32_t vu, int base_label,
                           ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources){
    // FUNC
    //     Given unitig vu, check if 
    //       1) it's a not-very-long tip that connects to exactly 2 vertices, 
    //       2) the 2 vertices is not a pair of haplotigs, and
    //       3) the 2 vertices could've been connected (i.e. the 2 ends connected with the tip overlap with each other)
    // RETURN
    //     1 if yes
    //     0 if no
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t wu[2];

    // tip direction: has no predecessor and has 2 successors
    if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)>0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=2){
        if (verbose){
            fprintf(stderr, "[debug::%s] utg %.6d , failed tip requirements\n", __func__, (int)(vu>>1)+1);
        }
        return 0;
    }

    // check haplotig status
    hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label);
    if (hamt_check_diploid(ug, wu[0], wu[1], 0.3, reverse_sources)>0){
        if (verbose){
            fprintf(stderr, "[debug::%s] utg %.6d , failed haplotig requirement\n", __func__, (int)(vu>>1)+1);
        }
        return 0;
    }

    // check overlap
    if (hamt_check_diploid(ug, wu[0], wu[1], 0, sources)>0){  // abuse the function: check homo paf with ratio threshold set to 0
        if (verbose){
            fprintf(stderr, "[debug::%s] utg %.6d , pass\n", __func__, (int)(vu>>1)+1);
        }
        return 1;
    }else{
        if (verbose){
            fprintf(stderr, "[debug::%s] utg %.6d , failed ovlp check\n", __func__, (int)(vu>>1)+1);
        }
        return 0;
    }
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
    asg_arc_t t_forward, t_rev, *p;
    int ql = coverage_cut[Get_qn(*h)].e - coverage_cut[Get_qn(*h)].s;
    int tl = coverage_cut[Get_tn(*h)].e - coverage_cut[Get_tn(*h)].s;
    int r_forward = ma_hit2arc(h, ql, tl, asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t_forward);
    ql = coverage_cut[Get_qn(*h_rev)].e - coverage_cut[Get_qn(*h_rev)].s;
    tl = coverage_cut[Get_tn(*h_rev)].e - coverage_cut[Get_tn(*h_rev)].s;
    int r_rev = ma_hit2arc(h_rev, ql, tl, asm_opt.max_hang_Len, asm_opt.max_hang_rate, asm_opt.min_overlap_Len, &t_rev);
    if (r_forward>=0 && r_rev>=0){
        if (yes_recover_it){ // sancheck: do nothing if arc exists
            fprintf(stderr, "[debug::%s]recover, forward %d , rev %d\n", __func__,  r_forward, r_rev);  // TODO/BUG: why r_forward!=r_rev ?
            p = asg_arc_pushp(sg);
            *p = t_forward;
            p = asg_arc_pushp(sg);
            *p = t_rev;
            if (verbose) fprintf(stderr, "[debug::%s] success\n", __func__);
            
        }
        return 1;
    }else{
        if (yes_recover_it && verbose) 
            fprintf(stderr, "[debug::%s] wanted to recover arc but gave up.\n", __func__);
        return 0;
    }    
}


int hamt_ug_recover_ovlp_if_existed(asg_t *sg, ma_ug_t *ug, uint32_t start, uint32_t end,
                                    ma_hit_t_alloc *sources, 
                                    const ma_sub_t* coverage_cut, int search_span,
                                    int yes_recover_it){
    // FUNC
    //     check if start->end could've formed an arc. Add the arc accordingly if so, or not (*set yes_recover_it)
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
    int verbose = 0;

    uint32_t v, w, v_start, v_end;
    asg_t *auxsg = ug->g;
    int status;

    v = start;
    w = end;
    status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, yes_recover_it);
    if (status<0){
        return -1;
    }else if (status>0){
        // fprintf(stderr, "ideal recovery\n");
        return 1;
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
            status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w_end, sources, coverage_cut, yes_recover_it);
            assert(status==1);
            // fprintf(stderr, "stepped recovery, read pair: %.*s - %.*s\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1),
            //                                                                (int)Get_NAME_LENGTH(R_INF, w_end>>1), Get_NAME(R_INF, w_end>>1));
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
        status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v, w, sources, coverage_cut, yes_recover_it);
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
            status = hamt_ug_recover_ovlp_if_existed_core(sg, ug, v_end, w, sources, coverage_cut, yes_recover_it);
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
int hamt_ug_check_complexBubble(asg_t *sg, ma_ug_t *ug, int max_size, uint32_t v0, uint32_t *end0, int base_label,
                                int max_length){
    // FUNC
    //    use a BFS variant to test if v0 is the start of a complex bubble (of size less than `max_size`)
    //    (complex bubble means: no inversion aka all vertices shall be discovered in only 1 direction,
    //                          no loop/circle)
    //    if so, arbitrarily chose a path from it and drop all the other arcs.
    //    If max_length>0, this function returns 0 if any of the bubble edge is long than max_length
    //    Another hidden criteria is any bubble edge shall not be significantly longer than the handle.
    // RETURN
    //    1 if yes (side effect: store the end vertex in *end0)
    //    0 if no
    // NOTE / TODO
    //    This implementation misses a few cases (can be easily handled by other routines), 
    //     but it can tolerates tips. 

    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, v_tmp, nv, nw, tmp_nv;
    asg_arc_t *av, *aw, *tmp_av;
    int nb_nodes_visited = 0;
    vu = v0;

    if (verbose) {fprintf(stderr, "[debug::%s] check vu utg%.6d\n", __func__, (int)(v0>>1)+1);}

    // basic topo checks, requiring vu to be a articulation-vertex-like vertex
    // (note that to this point, we should have no trivial tips.)
    // (forward direction)
    // NOTE: 
    //     disabling this check will select a single path many-tips-joined-together cases
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
                if ((tmp_av[i_tmp].v^1) == vu){continue;}  // note: must check direction
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

    int is_abnormal = 0, is_complex_bubble = 0, nb_rounds = 0;
    int has_anything_happend = 1;
    int v0_special_flag = 1;
    

    // test if we have a complex bubble here
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
                is_complex_bubble = 1;
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
    if (is_complex_bubble && (v0==(*end0))){
        if (ug->u.a[v0>>1].len<(longest_length_in_bubble*2)){
            is_complex_bubble = 0;
            if (verbose){
                fprintf(stderr, "[debug::%s] detected complex bubble for handle utg%.6d but discarded bc length (v0 length %d, longest %d)\n", __func__, (int)(v0>>1)+1, (int)ug->u.a[v0>>1].len, longest_length_in_bubble);
            }
        }
    }
    if (max_length>0 && longest_length_in_bubble>max_length){is_complex_bubble = 0;}

    ////////// debug sanchecks /////////////
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){  // color shouldn't be larger than 2
        if (color[i]>2){
            fprintf(stderr, "[debugERROR::%s] color can't be larger than 2 (utg%.6d)\n", __func__, (int)(i>>1)+1);fflush(stderr);
            return 0;// exit(1);
        }
    }
    ////////////////////////////////////////

    queue32_destroy(&q[0]);
    queue32_destroy(&q[1]);
    stacku32_destroy(&stack);
    free(color);

    return is_complex_bubble;
}


int hamt_ug_pop_subgraph(asg_t *sg, ma_ug_t *ug, uint32_t start0, uint32_t end0, 
                              int base_label, int alt_label){
    // TODO/BUG
    //     Has some issues in tangles; currently only used in complex bubble popping.
    // FUNC
    //     start0 and end0 is the start/end handle of a subgraph
    //     Will retain an arbitrary path that has max total reads
    //       and discard other arcs
    // NOTE
    //     Assumes that coverages have been collected.
    // RETURN
    //     number of arcs dropped
    int verbose = 0;
    asg_t *auxsg = ug->g;
    int nb_cut = 0;

    // (for pushing unitigs to alt hap)
    stacku32_t alt_buf;
    stacku32_init(&alt_buf);

    // get list of vertices involved in this complex bubble and do toposort
    stacku32_t vertices, finishing_times, DFSstack;
    stacku32_init(&DFSstack);
    stacku32_init(&vertices);
    stacku32_init(&finishing_times);
    uint64_t *packed;
    //  (get all vertices)
    stacku32_push(&DFSstack, start0);
    stacku32_push(&vertices, start0);
    stacku32_push(&finishing_times, 0);  // dummy
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    uint32_t vu;
    int is_exhausted;
    uint32_t nv;
    asg_arc_t *av;
    if (verbose>1){fprintf(stderr, "[debug::%s] > get the subgraph (%.6d - %.6d )\n", __func__, (int)(start0>>1)+1, (int)(end0>>1)+1);}
    while (stacku32_pop(&DFSstack, &vu)){
        is_exhausted = 1;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){ continue ;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if (av[i].v==end0){continue;}

            if (color[av[i].v]==0){
                is_exhausted = 0;
                color[av[i].v] = 1;
                stacku32_push(&DFSstack, vu);
                stacku32_push(&DFSstack, av[i].v);
                stacku32_push(&vertices, av[i].v);
                if (verbose>1){fprintf(stderr, "[debug::%s]     push %.6d\n", __func__, (int)(av[i].v>>1)+1);}
                stacku32_push(&finishing_times, 0);  // dummy
                break;
            }
        }
        if (is_exhausted){
            color[vu] = 2;
            if (verbose>1){fprintf(stderr, "[debug::%s]     ( finished %.6d )\n", __func__, (int)(vu>>1)+1);}
        }
    }
    // (push the end node)
    if (end0!=start0){
        color[end0] = 2;
        stacku32_push(&DFSstack, end0);
        stacku32_push(&vertices, end0);
        stacku32_push(&finishing_times, 0);  // dummy
    }

    //   (toposort)
    if (verbose){
        fprintf(stderr, "[debug::%s] start %.6d end %.6d size of %d\n", __func__,
                            (int)(start0>>1)+1, (int)(end0>>1)+1, vertices.n);
    }
    packed = (uint64_t*)malloc(sizeof(uint64_t) * vertices.n);
    asg_get_subgraph_DFSfinishTimes(auxsg, &vertices, &finishing_times);
    if (vertices.n<=1) {
        stacku32_destroy(&vertices);
        stacku32_destroy(&finishing_times);
        free(packed);

        return 0;
    }
    // fprintf(stderr, "%d, %d\n", vertices.n, finishing_times.n); fflush(stderr);
    assert(vertices.n==finishing_times.n);
    for (int i=0; i<vertices.n; i++){
        if (verbose>1){fprintf(stderr, "[debug::%s] %.6d finish t: %d\n", __func__, (int)(vertices.a[i]>>1)+1, (int)(finishing_times.a[i]));}
        packed[i] = ((uint64_t)finishing_times.a[i])<<32 | (uint64_t)vertices.a[i];
    }
    radix_sort_ovhamt64(packed, packed+vertices.n);
    //  (reverse the ordder: toposort order is the reverse to finishing time)
    uint64_t tmp;
    for (int i=0; i<vertices.n/2; i++){
        tmp = packed[i];
        packed[i] = packed[vertices.n-1-i];
        packed[vertices.n-1-i] = tmp;

    }
       
    // get the max coverage present
    int max_cov = 0;
    for (int i=0; i<vertices.n; i++){
        max_cov = ug->utg_coverage[vertices.a[i]>>1]>max_cov? ug->utg_coverage[vertices.a[i]>>1] : max_cov ;
    }
    max_cov+=1;
    if (verbose) fprintf(stderr, "[debug::%s]  (max coverage is %d)\n", __func__, max_cov);
    
    // find the path via shortest path algo
    //  (init)
    uint32_t *scores = (uint32_t*)malloc(sizeof(uint32_t) * vertices.n);
    uint32_t *pis = (uint32_t*)malloc(sizeof(uint32_t) * vertices.n);  // pi, aka the optimal path record
    for (int i=0; i<vertices.n; i++){
        scores[i] = 0; // 65535;
        pis[i] = 0;
    }

    // (traverse)
    uint32_t t, wu;
    int edgeweight;
    int i_wu=0, i_vu=0;
    //    (init)
    uint32_t v_start = (uint32_t)packed[0];
    for (int ii=0; ii<vertices.n; ii++){
        if (vertices.a[ii]==v_start) {
            i_vu = ii;
            break;
        }
    }
    assert(vertices.a[i_vu]==v_start);
    scores[i_vu] = 0;
    for (int i=0; i<finishing_times.n; i++){  // for each vertex in the sorted subgraph
        t = (uint32_t)(packed[i]>>32);
        vu = (uint32_t)packed[i];
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        if (verbose>1){fprintf(stderr, "[debug::%s] at i %d, utg %.6d\n", __func__, i, (int)(vu>>1)+1);}
        for (uint32_t j=0; j<nv; j++){  // for each target of the vertex
            wu = av[j].v;
            if (!stack32_is_in_stack(&vertices, wu)) continue;

            for (int ii=0; ii<vertices.n; ii++){
                if (vertices.a[ii]==wu) i_wu = ii;
                if (vertices.a[ii]==vu) i_vu = ii;
            }
            assert(vertices.a[i_wu]==wu);
            assert(vertices.a[i_vu]==vu);

            // edgeweight = max_cov - ug->utg_coverage[wu>>1];  // edge weight of vu=>wu; negating because we want the "longest" path
            edgeweight = ug->utg_coverage[wu>>1];
            if (verbose>1){fprintf(stderr, "[debug::%s]     target utg %.6d with edgecost %d (current score %d parent %.6d)\n", __func__, (int)(wu>>1)+1, edgeweight, scores[i_wu], (int)(pis[i_wu]>>1)+1);}

            // relax
            if (scores[i_wu]<scores[i_vu]+edgeweight){
                scores[i_wu] = scores[i_vu]+edgeweight;
                if (verbose>1){
                    fprintf(stderr, "[debug::%s]     ~ update linkage: %.6d 's parent changed to %.6d (score %d)\n", __func__, 
                    (int)(wu>>1)+1, (int)(vu>>1)+1, scores[i_wu]);
                }
                pis[i_wu] = vu;  // wu's "best" parent is vu
            }
        }
    }

    // pop
    //  (drop all arcs)
    // nb_cut = hamt_ug_arc_del_between(sg, ug, start0, end0, 1, base_label);
    //  (recover the selected path)
    int j=0, verbose_count=0;
    wu = end0;
    int encountered_error = 0;
    uint32_t nt, tu;
    asg_arc_t *at; 
    while (1){
        for (int i=0; i<vertices.n; i++){
            if (vertices.a[i]==wu) {
                i_wu = i;
                break;
            }
        }
        assert(vertices.a[i_wu]==wu);
        vu = pis[i_wu];  // arc is vu=>wu
        if (verbose>1){fprintf(stderr, "[debug::%s]  retaining %.6d -> %.6d \n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);}

        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].v!=wu) {
                hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
                stacku32_push(&alt_buf, av[i].v);
                nb_cut+=1;
            }else{
                tu = wu^1;
                nt = asg_arc_n(auxsg, tu);
                at = asg_arc_a(auxsg, tu);
                for (uint32_t j=0; j<nt; j++){
                    if (at[j].v!=(vu^1)){
                        hamt_ug_arc_del(sg, ug, tu, at[j].v, 1);
                        nb_cut+=1;
                        stacku32_push(&alt_buf, at[j].v);
                    }
                }
            }
        }

        if (vu==start0) break;

        wu = vu;
        verbose_count+=1;
        if (verbose_count>1000) break;  // DEBUG

    }
    while (stacku32_pop(&alt_buf, &vu)){
        hamt_ug_utg_softdel(sg, ug, vu, alt_label);
    }

    if (verbose){fprintf(stderr, "[debug::%s] retained path has %d arcs\n", __func__, verbose_count);}


    stacku32_destroy(&vertices);
    stacku32_destroy(&finishing_times);
    stacku32_destroy(&alt_buf);
    free(packed);
    free(scores); free(pis);

    return nb_cut;
}

int hamt_ug_pop_complexBubble(asg_t *sg, ma_ug_t *ug, uint32_t start0, uint32_t end0, 
                              int base_label, int alt_label, int is_hard_drop){
    // FUNC
    //     start0 and end0 is the start/end handle of a complex bubble (per `hamt_ug_check_complexBubble`)
    //     retain an arbitrary path (greedily long but not the longest) and discard other arcs
    // TODO
    //     try to be not too arbitrary
    // RETURN
    //     number of dropped arcs
    int verbose = 0;

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
int hamt_ug_pop_complexBubble_v2(asg_t *sg, ma_ug_t *ug, uint32_t start0, uint32_t end0, 
                              int base_label, int alt_label, int is_hard_drop){
    // RETURN
    //     number of arcs dropped
    int ret = hamt_ug_pop_subgraph(sg, ug, start0, end0, base_label, alt_label);
    return ret;
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
    // hamt_asgarc_ugTreatBifurcation(sg, ug, -1, 10000000, base_label);
}

int hamt_ug_pop_bubble(asg_t *sg,ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop)
{
    // remove simple bubbles at unitig graph level
    int verbose = 0;
    double startTime = Get_T();
    int max_edge_length = 300000;  // don't touch large bubble-like structure on unitig graph

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
        
        if (ug->u.a[wu[0]>>1].len>max_edge_length || ug->u.a[wu[1]>>1].len>max_edge_length){
            continue;
        }

        // // pop (length or coverage)
        // if (ug->u.a[wu[0]>>1].len>50000 || ug->u.a[wu[1]>>1].len>50000){
        //     // if either side of the bubble is long,
        //     // pop the shorter side
        //     idx = ug->u.a[wu[0]>>1].len > ug->u.a[wu[1]>>1].len ? 1 : 0;
        // }else{
        //     // otherwise, pop the side with lesser coverage
        //     idx = ug->utg_coverage[wu[0]>>1] > ug->utg_coverage[wu[1]>>1]? 1:0;
        // }

        // pop the side with lesser coverage
        idx = ug->utg_coverage[wu[0]>>1] > ug->utg_coverage[wu[1]>>1]? 1:0;  // idx is the one to be dropped
        
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
        if (idx<2 || idx>=20){
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
    }
    return nb_cut;
}

int hamt_ug_drop_transitive_links(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //     Given a->b->c and a->c (and the reverse direction) and a!=c, 
    //      we consider drop a->c (and the rev).
    //     (in the code a,b,c is vu,wu,uu)
    // TODO
    //     Only checking the most obvious cases right now (you can only have one b).
    //     The ideal way is to DFS with limits of how many nodes/base pairs for the shortest path.
    // RETURN
    //     Number of arcs dropped.
    int ret = 0;
    int verbose = 0;
    uint32_t wu, uu, nv, nu;
    int yes_del;
    asg_arc_t *av, *au;
    asg_t *auxsg = ug->g;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (auxsg->seq_vis[vu>>1]!=base_label) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<=1) continue;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            wu = av[i].v;
            if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)<=1) continue;

            yes_del = 0;  // reset flag
            for (int j=0; j<nv; j++){
                if (j==i) continue;  // it's wu it self
                if (av[j].del) continue;
                uu = av[j].v;
                if ((uu>>1)==(vu>>1)) continue;
                nu = asg_arc_n(auxsg, uu);
                au = asg_arc_a(auxsg, uu);
                for (int z=0; z<nu;z++){
                    if (au[z].del) continue;
                    if (au[z].v==wu){
                        if (verbose){
                            fprintf(stderr, "[debug::%s] drop arc between utg%.6d and utg%.6d\n",
                                            __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
                        }
                        yes_del = 1;
                        ret++;
                        hamt_ug_arc_del(sg, ug, vu, wu, 1);
                        break;
                    }
                }
                if (yes_del) break;  // vu->wu has been removed, no need to check other paths.
            }
        }
    }
    if (ret){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return ret;
}
int hamt_ug_drop_transitive_nodes(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label){
    // FUNC
    //     Remove a node if all its predecessors could still easily reach all its
    //      targets after removal.
    // PARAMETER
    //     Do not remove nodes larger than `size_limit_bp`, since
    //      we're only checking the immediate surroundings and therefore
    //      not aware of loops.
    // TODO
    //     We are only checking the most obvious case right now, 
    //       i.e. the predecessors shall reach the targets via passes
    //       only containing one additional node.
    // RETURN
    //     Number of nodes dropped.
    int ret = 0;
    int verbose = 0;
    uint32_t vu_pre, vu_tar;
    int yes_del;

    uint32_t n_pre, n_tar;
    asg_arc_t *a_pre, *a_tar;
    asg_t *auxsg = ug->g;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>size_limit_bp) continue;
        if (auxsg->seq_vis[vu>>1]!=base_label) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==0 ||
            hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)==0) continue;  // tip or free floating unitig
        n_pre = asg_arc_n(auxsg, vu^1);
        a_pre = asg_arc_a(auxsg, vu^1);
        n_tar = asg_arc_n(auxsg, vu);
        a_tar = asg_arc_a(auxsg, vu);

        if (verbose){
            fprintf(stderr, "[debug::%s] at utg%.6d. Predecessors:\n", __func__, (int)(vu>>1)+1);
            for (uint32_t i_pre=0; i_pre<n_pre; i_pre++){
                if (a_pre[i_pre].del) continue;
                fprintf(stderr, "[debug::%s]     utg%.6d\n", __func__, (int)(a_pre[i_pre].v>>1)+1);
            }
            fprintf(stderr, "[debug::%s] (Successors:)\n", __func__);
            for (uint32_t i_tar=0; i_tar<n_tar; i_tar++){
                if (a_tar[i_tar].del) continue;
                fprintf(stderr, "[debug::%s]     utg%.6d\n", __func__, (int)(a_tar[i_tar].v>>1)+1);
            }
        }

        yes_del = 1;  // reset
        for (uint32_t i_pre=0; i_pre<n_pre; i_pre++){
            if (a_pre[i_pre].del) continue;
            vu_pre = a_pre[i_pre].v^1;
            for (uint32_t i_tar=0; i_tar<n_tar; i_tar++){
                if (a_tar[i_tar].del) continue;
                vu_tar = a_tar[i_tar].v;
                if (hamt_check_if_jump_reachable_simple(auxsg, vu_pre, vu_tar, vu, 1, 0)==0){
                    yes_del = 0;
                    break;
                }
            }
            if (!yes_del) break;  // at least one pair isn't reacheable 
        }
        if (yes_del){
            if (verbose){
                fprintf(stderr, "[debug::%s] drop utg%.6d\n", __func__, (int)(vu>>1)+1);
            }
            for (uint32_t i_pre=0; i_pre<n_pre; i_pre++){
                hamt_ug_arc_del(sg, ug, a_pre[i_pre].v^1, vu, 1);
            }
            for (uint32_t i_tar=0; i_tar<n_tar; i_tar++){
                hamt_ug_arc_del(sg, ug, vu, a_tar[i_tar].v, 1);
            }
            auxsg->seq_vis[vu>>1] = !base_label;
            ret++;
        }

    }
    if (ret){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return ret;
}
int hamt_ug_drop_transitive(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label){
    // A better hamt_ug_pop_miscbubble_aggressive.
    // RETURN
    //     Total number of treatments. 

    // Drop links first. 
    // The other way around is ok too, but note that the result will be slightly differenct.
    int ret = 0;
    for (int i=0; i<3; i++){
        int ret1=0, ret2=0;
        ret1 = hamt_ug_drop_transitive_links(sg, ug, base_label);

        // Drop nodes. 
        ret2 = hamt_ug_drop_transitive_nodes(sg, ug, size_limit_bp, base_label);

        fprintf(stderr, "[M::%s] dropped %d links and %d nodes.\n", __func__, ret1, ret2);
        ret += ret1+ret2;
        
        if (ret1+ret2==2)break;
    }

    return ret;
}


int hamt_ug_pop_simpleInvertBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop){
    int verbose = 0;
    double startTime = Get_T();

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu;
    int nb_cut = 0;

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>100000){
            continue;
        }
        if (verbose>1){fprintf(stderr, "[debug::%s] at utg%.6d\n", __func__, (int)(vu>>1)+1);}
        if (hamt_asgarc_util_checkSimpleInvertBubble(auxsg, vu, &wu, &uu, base_label)){  // note: uu's direction is adjusted for popping; no need to ^1
            if (ug->u.a[wu>>1].len>100000){
                continue;
            }
            // if (ug->u.a[uu>>1].len>50000){
            //     continue;
            // }
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
                fprintf(stderr, "[debug::%s] pop, v utg%.6d w utg%.6d u utg%.6d \n", __func__, 
                                (int)(vu>>1)+1, (int)(wu>>1)+1, (int)(uu>>1)+1);
            }
        }
    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    // if (VERBOSE){
        fprintf(stderr, "[M::%s] popped %d locations\n", __func__, nb_cut);
    // }
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

int hamt_ug_pop_tinyFlatCircles(asg_t *sg, ma_ug_t *ug, int base_label){
    /*
    ------v0   v2-------
           \  /
            v1(a self circle unitig)
     (v0 could equals v2)
     (v0 and v2 can have other linkage)
     If v1's coverage is similar to v0 and v2, then drop the 
       self-circle arc for v1.
    */
   // RETURN
   //     count of arcs dropped.
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t v0, v1, v2, vtmp[2], w1, w2;
    int cov0, cov1, cov2, diffa, diffb;
    int total = 0;

    for (v1=0; v1<auxsg->n_seq; v1++){
        if (hamt_asgarc_util_countSuc(auxsg, v1, 0, 0, base_label)!=2 || 
            hamt_asgarc_util_countPre(auxsg, v1, 0, 0, base_label)!=2
            ){
                continue;
            }
        hamt_asgarc_util_get_the_two_targets(auxsg, v1, &vtmp[0], &vtmp[1], 0, 0, base_label);
        if (vtmp[0]!=v1 && vtmp[1]!=v1){continue;}  // v1 is not a self circle
        if (vtmp[0]==v1 and vtmp[1]==v1){  // sancheck, should not happen
            fprintf(stderr, "[E::%s] double arc at %.6d ? continue anyway\n", __func__, (int)(v1>>1));
            continue;
        }
        v0 = vtmp[0]==v1? vtmp[1] : vtmp[0];

        hamt_asgarc_util_get_the_two_targets(auxsg, v1^1, &vtmp[0], &vtmp[1], 0, 0, base_label);
        if (vtmp[0]!=(v1^1) && vtmp[1]!=(v1^1)){ // should not happen, v1 is not a normal self circle
            fprintf(stderr, "[E::%s] asym arc? %.6d\n", __func__, (int)(v1>>1));
            continue;
        }  
        if (vtmp[0]==(v1^1) and vtmp[1]==(v1^1)){  // sancheck, should not happen
            fprintf(stderr, "[E::%s] double arc at %.6d ? continue anyway\n", __func__, (int)(v1>>1));
            continue;
        }
        v2 = vtmp[0]==(v1^1)? vtmp[1] : vtmp[0];

        // collect coverage
        cov0 = ug->utg_coverage[v0>>1];
        cov1 = ug->utg_coverage[v1>>1];
        cov2 = ug->utg_coverage[v2>>1];
        diffa = cov1>cov0? cov1-cov0 : cov0-cov1;
        w1 = cov1>cov0? cov1 : cov0;
        diffb = cov2>cov0? cov2-cov0 : cov0-cov2;
        w2 = cov2>cov0? cov2 : cov0;

        // check coverage diff
        if ((float)diffa/w1 > 0.3 || (float)diffb/w2 > 0.3){
            continue;
        }

        // cut
        total++;
        if (verbose>1){
            fprintf(stderr, "[debug::%s] v0 %.6d (%d)  v1 %.6d (%d)  v2 %.6d (%d)\n", __func__,
                                (int)(v0>>1)+1, cov0,
                                (int)(v1>>1)+1, cov1,
                                (int)(v2>>1)+1, cov2);
        }
        hamt_ug_arc_del_selfcircle(sg, ug, v1, base_label);
        if ((v1>>1)==16){  // debug
            uint32_t nv;
            asg_arc_t *av;
            nv = asg_arc_n(auxsg, v1);
            av = asg_arc_a(auxsg, v1);
            for (int i=0; i<nv; i++){
                fprintf(stderr, "    (1) target %.6d del %d\n", (int)(av[i].v>>1)+1, (int)av[i].del);
            }
        }
    }

    if (verbose){
        fprintf(stderr, "[debug::%s] treated %d spots\n", __func__, total);
    }
    if (total){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return total;
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
        idx = -1;
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
        if (idx<0) {continue;}  // shouldn't happen
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
    // by one or multiple arc(s) shooting from ONE end /(the target also has only this one predecessor)/,
    // cut it off regardless of coverage diff
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t vu, vu_, wu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=1){continue;}
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)<2/*!=2*/){continue;}
        if (hamt_asgarc_util_get_the_one_target(auxsg, vu, &vu_, 0, 0, base_label)<0){fprintf(stderr, "ERROR %s\n", __func__); fflush(stderr); exit(0);}
        if (vu_!=vu){continue;}  // check circling
        nv = asg_arc_n(auxsg, vu^1);
        av = asg_arc_a(auxsg, vu^1);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if ((av[i].v^1)==vu){continue;}
            wu = av[i].v;
            if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, base_label)==0){  // spare if is long tip
                if (ug->u.a[wu>>1].len>100000){
                    continue;
                }
            }  
            // cut
            hamt_ug_arc_del(sg, ug, vu^1, wu, 1);  // compile note: safe by assertion
            nb_cut++;
            if (verbose){fprintf(stderr, "[M::%s] cut between utg%.6d and utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
        }
        }

    }

    fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
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
        Will cut the vu->wu arc. Will send vu to the alternative collection.
    */
    asg_t *auxsg = ug->g;
    uint32_t wu, uu=0, nv;
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
            hamt_ug_utg_softdel(sg, ug, vu, base_label^1);
            nb_cut++;
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    // if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d mid length tips.\n", __func__, nb_cut);
    // }
    return nb_cut;
}
int hamt_ug_drop_midsizeTips_aggressive(asg_t *sg, ma_ug_t *ug, float fold, int base_label){
    // FUNC
    //     treat the follow tip
    /*
               vu
              /
    --------wu----uu--------
        difference: 
            - vu is strictly a simple tig; wu and uu can have other connections
            - drop vu if vu is shorter than wu and uu
    */
    asg_t *auxsg = ug->g;
    uint32_t wu, uu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    float r = fold<0? 0.5 : fold;
    int length_check_passed = 0;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        length_check_passed = 0;
        // tip topo
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=1){continue;}

        // wu
        if (hamt_asgarc_util_get_the_one_target(auxsg, vu, &wu, 0, 0, base_label)<0){
            fprintf(stderr, "[E::%s] can't get target\n", __func__);
            continue;
        }
        if (hamt_asgarc_util_countPre(auxsg, wu, 0, 0, base_label)<2){continue;}


        if ( !(ug->u.a[wu>>1].len>(ug->u.a[vu>>1].len*r)) ){
            continue;
        }

        // uu
        nv = asg_arc_n(auxsg, wu^1);
        av = asg_arc_a(auxsg, wu^1);
        uint32_t i;
        for (i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (base_label>=0 && auxsg->seq_vis[av[i].v>>1]!=base_label) {continue;}
            if ((av[i].v>>1)==(vu>>1)){continue;}
            uu = av[i].v;
            if (ug->u.a[uu>>1].len>(ug->u.a[vu>>1].len*r)){
                length_check_passed = 1;
                break;
            }
        }

        if (length_check_passed){
            // cut
            hamt_ug_arc_del(sg, ug, vu, wu, 1);
            hamt_ug_utg_softdel(sg, ug, vu, base_label^1);
            nb_cut++;
        }
    }
    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    // if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d mid length tips.\n", __func__, nb_cut);
    // }
    return nb_cut;
}

void hamt_ug_prectgTopoClean(asg_t *sg, 
                            const ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex,
                            int base_label, int alt_label, int is_hard_drop){
    // note: will redo the unitig graph
    // assume graph is clean
    // pop all simple bubbles, apparent bubbles, pentagon bubbles (popping doesn't depend on coverage)
    //  (length requirement is arbitrary)
    double startTime = Get_T();
    int verbose = 0;

    ma_ug_t *ug = hamt_ug_gen(sg, coverage_cut, sources, ruIndex, base_label);
    asg_t *auxsg = ug->g;
    uint32_t vu, wu, uu, nv;
    asg_arc_t *av;
    int nb_pop = 0, round = 0;

    for (round=0; round<3; round++){
        if (asm_opt.write_debug_gfa) {hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, "prectg-TOPO", round);}

        nb_pop = 0;

        // inverted simple bubbles
        nb_pop += hamt_ug_pop_simpleInvertBubble(sg, ug, base_label, alt_label, is_hard_drop);

        // penta simple bubbles
        uint32_t penta_buf[3];
        int nb_penta_pop = 0;
        for (vu=0; vu<auxsg->n_seq*2; vu++){
            if (hamt_asgarc_util_checkSimplePentaBubble(ug, auxsg, vu, penta_buf, base_label)>0){

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
            fprintf(stderr, "[M::%s] popped total %d locations (round %d).\n", __func__, nb_pop, round);
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
    //    Only use after basic topo clean has been done, this function will only check very simple cases.
    //    Does several rounds of ug destroy-regen, also writes inter gfa if specified.
    int verbose = 0;
    int nb_modified = 0;
    double startTime;

    uint32_t start, end, start_v, end_v;  // start and end are meant to be unitig IDs (with dir); start_v and end_v are vertex IDs
    int l1, l2;
    uint32_t nv, wu, uu;
    asg_arc_t *av;
    int is_circular;

    for (int round=0; round<3; round++){
        startTime = Get_T();
        nb_modified = 0;
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
                if (verbose && l1>0){
                    fprintf(stderr, "[W::%s] something wrong 1a, vu is utg%.6d, start utg%.6d, end utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(start>>1)+1, (int)(end>>1)+1);
                }
                continue;
            }
            if (hamt_asgarc_util_countPre(auxsg, start, 0, 0, base_label)!=1){
                if (verbose && l2>0){
                    fprintf(stderr, "[W::%s] something wrong 1b, vu is utg%.6d, start utg%.6d, end utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(start>>1)+1, (int)(end>>1)+1);
                }
                continue;
            }

            // check if after stepping 1 step, the end unitg will find the start unitig
            if (verbose){
                fprintf(stderr, "[debug::%s]     got to check if there was any overlap\n", __func__);
            }
            if (hamt_asgarc_util_get_the_one_target(auxsg, end, &wu, 0, 0, base_label)<0){
                if  (verbose) {fprintf(stderr, "something wrong 2%s\n", __func__);}
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
                    if (hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, 5, 1)>0){
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
        if (VERBOSE){
            fprintf(stderr, "[M::%s] finished round %d, modified %d, used %.2f s\n", __func__, round, nb_modified, Get_T()-startTime);
        }
        // hamtdebug_output_unitig_graph_ug(ug, asm_opt.output_file_name, 210+round);
        if (nb_modified){
            free(sg->idx);
            sg->idx = 0;
            sg->is_srt = 0;
            asg_cleanup(sg);
            asg_cleanup(auxsg);
            hamt_ug_destroy(ug);
        }else{
            hamt_ug_destroy(ug);
            if (VERBOSE) {fprintf(stderr, "[M::%s] leaving\n", __func__);}
            break;
        }
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
                    if (hamt_ug_recover_ovlp_if_existed(sg, ug, w, v, sources, coverage_cut, 5, 1)>0){
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
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_treat);
    }

}

void hamt_ug_rescueLongUtg(asg_t *sg, 
                                    ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,R_to_U* ruIndex,
                                    const ma_sub_t* coverage_cut)
{
    // experimental
    //    if a long unitig is not connected and non-circular,
    //    try the ends and add a link if there was an overlap in the paf
    //    Obviously more aggressive than hamt_ug_prectg_rescueShortCircuit, which required more topo checks. 
    int verbose = 0;

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
        if (hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, 20, 1)>0){
            nb_treated++;
            color[vu>>1] = 1;
            if (verbose){
                fprintf(stderr, "[debug::%s] added circle link for utg%.6d (homo)\n", __func__, (int)(vu>>1)+1);
            }
        }else if (hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, reverse_sources, coverage_cut, 20, 1)>0){
            nb_treated++;
            color[vu>>1] = 1;
            if (verbose){
                fprintf(stderr, "[debug::%s] added circle link for utg%.6d (het)\n", __func__, (int)(vu>>1)+1);
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
int hamt_ug_try_circularize(asg_t *sg, ma_ug_t *ug, 
                            ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,R_to_U* ruIndex,
                            const ma_sub_t* coverage_cut, int l_threshold){
    // FUNC
    //     Check all long contigs with at least 1 arc for start-end overlap.
    //     If an overlap looks solid and the start & end is connected via other small contigs,
    //      consider to make it a circle and cut off from the current subgraph.
    //     For example:
    //            v1(start)---v2(start)---...(other contigs)
    //             |           |
    //            v1(end)-----v2(end)---...(other contigs)
    //        where v1 is long and v2 is tiny. There's no topo requirement for other arcs of v2.
    // NOTE
    //      To try to circularize isolated long linear contigs, see hamt_ug_rescueLongUtg.
    // RETURN
    //      Number of treatments
    int verbose = 0;
    int ret = 0;
    asg_t *auxsg = ug->g;
    uint32_t start_v, end_v;  // start/end read of unitig
    uint32_t start, target, tmp;
    uint32_t nv; 
    asg_arc_t *av;

    int search_span = 10;  // for finding overlapped reads at ends
    int search_neighbour_span = 3;  // for naive short path finding
    int search_neighbour_bp = 50000;  // for naive short path finding, max in-between contig length
    int search_neighbour_has_long = 0;  // for naive short path finding, if any connected contig is long
    queue32_t q[2];
    queue32_init(&q[0]);
    queue32_init(&q[1]);
    uint8_t which_q = 0;
    int passed;
    int sancheck_nbcut;
    
    int existed_homo = 0, existed_het = 0;
    int will_need_adding_selfarc = 0;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){  

        // candiate unitig should be long and not isolated
        if (ug->u.a[vu>>1].len<1000000) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0 && hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0) continue;
        
        // ignore unconnected circular contig; check if self arc already exists to avoid double arc
        will_need_adding_selfarc = hamt_asgarc_util_has_self_arc(auxsg, vu, 0);
        if (will_need_adding_selfarc==1) continue;  // unconnected circular contig
        will_need_adding_selfarc = !will_need_adding_selfarc;

        start_v = ug->u.a[vu>>1].start;
        end_v = ug->u.a[vu>>1].end^1;

        passed = 0;
        // check for overlap between start/end of the contig, and the shortest path connecting them
        existed_homo = hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, search_span, 0);
        existed_het = hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, reverse_sources, coverage_cut, search_span, 0);
        if (existed_homo>0 || existed_het>0){  // only check, don't push anything
            if (verbose){ fprintf(stderr, "[debug::%s] tig %.6d start %.*s end %.*s\n", __func__,
                                    (int)(vu>>1)+1, 
                                    (int)Get_NAME_LENGTH(R_INF, start_v>>1), 
                                         Get_NAME(R_INF, start_v>>1),
                                    (int)Get_NAME_LENGTH(R_INF, end_v>>1), 
                                         Get_NAME(R_INF, end_v>>1));
                        }
            // check if there's a very short path 
            // (naive, short in terms of both base pairs and number of contigs)
            if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0){
                target = vu^1;
                start = vu^1;
            }else{
                target = vu;
                start = vu;
            }
            // (init)
            queue32_reset(&q[0]);
            queue32_reset(&q[1]);
            search_neighbour_has_long = 0;
            nv = asg_arc_n(auxsg, start);
            av = asg_arc_a(auxsg, start);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if (av[i].v==target){
                    passed = 2;
                    break;
                }
                if (ug->u.a[av[i].v>>1].len>search_neighbour_bp){
                    search_neighbour_has_long = 1;  // don't push long neighbours
                }else{
                    queue32_enqueue(&q[which_q], av[i].v);
                }
            }
            // if no immediate arc, check a couple neighbours
            if (!passed){ 
                for (int i=0; i<search_neighbour_span; i++){
                    while(queue32_dequeue(&q[which_q], &tmp)){
                        if (tmp==target){
                            passed = 2;
                            break;
                        }
                        nv = asg_arc_n(auxsg, start);
                        av = asg_arc_a(auxsg, start);
                        for (int tmpi=0; tmpi<nv; tmpi++){  // push child nodes
                            if (av[i].del) continue;
                            if (ug->u.a[av[i].v>>1].len>search_neighbour_bp){
                                search_neighbour_has_long = 1;  // don't push long neighbours
                            }else{
                                queue32_enqueue(&q[!which_q], av[i].v);  // put neighbors in the other buffer
                            }
                        }
                    }
                    if (passed>1){break;}
                    queue32_reset(&q[which_q]);
                    which_q = !which_q;
                }
            }
            // do we want to circularize?
            if (verbose) {fprintf(stderr, "[deubg::%s]   pass %d, long neighbor %d\n", __func__, passed, search_neighbour_has_long);}
            if (passed>1 || (!search_neighbour_has_long)){
                sancheck_nbcut = 0;
                // drop links other than self link (if it exists)
                nv = asg_arc_n(auxsg, vu);
                av = asg_arc_a(auxsg, vu);
                for (int i=0; i<nv; i++){
                    if (av[i].v!=vu) {
                        hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
                        sancheck_nbcut++;
                    }
                }
                nv = asg_arc_n(auxsg, vu^1);
                av = asg_arc_a(auxsg, vu^1);
                for (int i=0; i<nv; i++){
                    if (av[i].v!=(vu^1)) {
                        hamt_ug_arc_del(sg, ug, vu^1, av[i].v, 1);
                        sancheck_nbcut++;
                    }
                }
                if (verbose) {fprintf(stderr, "[deubg::%s]   dropped %d, will=%d\n", __func__, sancheck_nbcut, will_need_adding_selfarc);}
                if (existed_homo>0){
                    hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, sources, coverage_cut, search_span, will_need_adding_selfarc);
                }else{
                    hamt_ug_recover_ovlp_if_existed(sg, ug, end_v, start_v, reverse_sources, coverage_cut, search_span, will_need_adding_selfarc);
                }
                ret++;
            }
        }
    }

    queue32_destroy(&q[0]);
    queue32_destroy(&q[1]);
    if (ret){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] added %d links.\n", __func__, ret);
    }
    return ret;
}

int hamt_ug_prectg_resolve_complex_bubble(asg_t *sg, ma_ug_t *ug, 
                                          int base_label, int alt_label, int is_hard_drop,
                                          int max_length){
    int verbose = 0;

    // ma_ug_t *ug = ma_ug_gen(sg);
    asg_t *auxsg = ug->g;
    uint32_t vu, nv, end, wu;
    asg_arc_t *av;
    int nb_cut = 0, nb_bubbles = 0;

    vecu64_t vec, vec_seen;
    vecu64_init(&vec);
    vecu64_init(&vec_seen);

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        // (expects vu to be the start of a complex bubble)
        if (hamt_ug_check_complexBubble(sg, ug, 500, vu, &end, base_label, max_length)>0){
            vecu64_push(&vec, ( ((uint64_t)vu)<<32 )|( (uint64_t)end ) );
            nb_bubbles++;
        }
    }

    uint64_t key, key_rev;
    if (nb_bubbles){
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
            nb_cut += hamt_ug_pop_complexBubble_v2(sg, ug, vu, end, base_label, alt_label, is_hard_drop);   
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
        if (tip_round>10){break;}  // waterproof just in case
    }
    asg_cleanup(sg);
    asg_cleanup(auxsg);

    fprintf(stderr, "[M::%s] dropped %d arcs for %d complex bubbles\n", __func__, nb_cut, nb_bubbles);
    return nb_cut;
}

int hamt_ug_resolve_oneMultiLeafSoapBubble(asg_t *sg, ma_ug_t *ug, 
                                           int base_label, int alt_label, int is_hard_drop){
    int nb_cut = 0;
    for (uint32_t vu=0; vu<ug->g->n_seq*2; vu++){
        nb_cut += hamt_ug_checknpop_oneMultiLeafSoapBubble(sg, ug, vu, base_label, alt_label, is_hard_drop);
    }
    // if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_cut);
    // }
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
    // if (VERBOSE){
        fprintf(stderr, "[M::%s] dropped %d arcs\n", __func__, nb_cut);
    // }
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
    nb_cut += hamt_ug_pop_miscbubble(sg, ug, base_label);
    nb_cut += hamt_ug_pop_simpleInvertBubble(sg, ug, base_label, alt_label, is_hard_drop);
    // nb_cut += hamt_ug_pop_miscbubble_aggressive(sg, ug, base_label);
    // nb_cut += hamt_ug_drop_transitive(sg, ug, 100000, base_label);

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
    fprintf(stderr, "[M::%s] total cut: %d\n", __func__, nb_cut);
    return nb_cut;
}

int hamt_ug_special_toploclean(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label){
    int nb_cut = 0;
    nb_cut += hamt_ug_popSpecialCase1(sg, ug, base_label, alt_label);
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
    int verbose = 0;

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
                    ret = hamt_ug_recover_ovlp_if_existed(sg, ug, start, end, reverse_sources, coverage_cut, 10, 1);
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
                ret = hamt_ug_recover_ovlp_if_existed(sg, ug, start, end, reverse_sources, coverage_cut, 10, 1);
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
    int verbose = 0;
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
    int verbose = 0;

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

void hamt_debug_dump(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources){
    asg_t *auxsg = ug->g;

    // a quick check of overlaps between the two ends of a given long linear unitig,
    //  which only examines the start and the end vertices. (modify hamt_ug_recover_ovlp_if_existed if need to search a span)
    // (note: any long linear unitig, not neceessarily unconnected)
    uint32_t start, end;
    ma_hit_t *ph1, *ph2;  // placeholder
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1){continue;}  // only need to check in one direction
        if (ug->u.a[vu>>1].len<500000){continue;}
        if (ug->u.a[vu>>1].circ){continue;}
        start = ug->u.a[vu>>1].start;
        end = ug->u.a[vu>>1].end^1;
        if (does_ovlp_ever_exist(sources, start, end, &ph1, &ph2)){
            fprintf(stderr, "[dbg-linearSancheck] tg%.6d cis ovlp\n", (int)(vu>>1)+1);
        }else if  (does_ovlp_ever_exist(reverse_sources, start, end, &ph1, &ph2)){
            fprintf(stderr, "[dbg-linearSancheck] tg%.6d trans ovlp\n", (int)(vu>>1)+1);
        }else{
            fprintf(stderr, "[dbg-linearSancheck] tg%.6d no ovlp\n", (int)(vu>>1)+1);
        }
    }
    fprintf(stderr, "\n\n");

    // haplotig info
    //    (only check the given tig's 2-level neighbours (undirected))
    vecu32_t buf;
    vecu32_init(&buf);
    uint32_t nv, nw, wu;
    asg_arc_t *av, *aw;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1){continue;}  // only need to check in one direction
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0){
            continue;
        }
        if (ug->u.a[vu>>1].circ){continue;}
        vecu32_reset(&buf);
        vecu32_push(&buf, vu);
        vecu32_push(&buf, vu^1);

        // 1st direction
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (auxsg->seq_vis[vu>>1]!=0){continue;}
            wu = av[i].v;
            if (vecu32_is_in_vec(&buf, wu)){continue;}
            vecu32_push(&buf, wu);
            // push wu's children
            nw = asg_arc_n(auxsg, wu);
            aw = asg_arc_a(auxsg, wu);
            for (uint32_t j=0; j<nw; j++){
                if (aw[j].del) {continue;}
                if (auxsg->seq_vis[wu>>1]!=0){continue;}
                wu = aw[j].v;
                if (vecu32_is_in_vec(&buf, wu)){continue;}
                vecu32_push(&buf, wu);
            }
        }
        // 2nd direction
        nv = asg_arc_n(auxsg, vu^1);
        av = asg_arc_a(auxsg, vu^1);
        for (uint32_t i=0; i<nv; i++){  // (same as above)
            if (av[i].del) {continue;}
            if (auxsg->seq_vis[vu>>1]!=0){continue;}
            wu = av[i].v;
            if (vecu32_is_in_vec(&buf, wu)){continue;}
            vecu32_push(&buf, wu);
            // push wu's children
            nw = asg_arc_n(auxsg, wu);
            aw = asg_arc_a(auxsg, wu);
            for (uint32_t j=0; j<nw; j++){
                if (aw[j].del) {continue;}
                if (auxsg->seq_vis[wu>>1]!=0){continue;}
                wu = aw[j].v;
                if (vecu32_is_in_vec(&buf, wu)){continue;}
                vecu32_push(&buf, wu);
            }
        }

        // check haplotig status between vu and these nodes
        fprintf(stderr, "[dbg-checkDip] base tg%.6d (dir %d)\n", (int)(vu>>1)+1, (int)(vu&1));
        for (int i=0; i<buf.n; i++){
            if ((buf.a[i]>>1)==(vu>>1)){continue;}
            fprintf(stderr, "[dbg-checkDip]   target tg%.6d (dir %d), ", (int)(buf.a[i]>>1)+1, (int)(buf.a[i]&1));
            if (hamt_check_diploid(ug, vu, buf.a[i], 0.5, reverse_sources)){  // note: ratio applies to the shorter tig
                fprintf(stderr, "yes\n");
            }else{
                fprintf(stderr, "no\n");
            }
        }

    }
    vecu32_destroy(&buf);

}

int hamt_ug_popTangles(asg_t *sg, ma_ug_t *ug, uint32_t source, uint32_t sink, 
                       uint32_t *buf_blacklist,  // note: maybe we don't need this. test more
                       int base_label, int alt_label){
    // NOTE
    //     direction is source->...<-sink
    // TODO
    //     Greedily try to find the path with more reads by
    //       searching the vertices with more coverages first.
    // RETURN
    //     1 if popped something
    //     0 otherwise
    int verbose = 0;
    uint32_t vu, nv, wu, uu, nu;
    asg_arc_t *av;
    asg_arc_t *au;
    asg_t *auxsg = ug->g;
    int ret = 0;

    // special case, early termination: self circle (shortcut or 1 unitig circle)
    nv = asg_arc_n(auxsg, source);
    av = asg_arc_a(auxsg, source);
    for (int i=0; i<nv; i++){
        if (av[i].del) continue;
        if (av[i].v==source){
            if (verbose) {fprintf(stderr, "[debug::%s] special case failure\n", __func__);}
            return 0;
        }
    }

    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq, 1);

    // a buffer for greedily choosing the path with most coverage
    stacku64_t for_sorting;
    stacku64_init(&for_sorting);
    uint64_t packed_value;

    stacku32_t s, visited;
    stacku32_init(&s);
    stacku32_init(&visited);  // directed IDs; linear search
    int is_exhausted, is_loop, san=0;
    int is_hinted_circular = ((sink>>1) == (source>>1));
    int is_reached_sink = 0;

    stacku32_push(&s, source);
    stacku32_push(&visited, source);
    stacku32_push(&visited, source^1);
    stacku32_push(&visited, sink);
    stacku32_push(&visited, sink^1);
    while (stacku32_pop(&s, &vu)){
        if (verbose) {fprintf(stderr, "[debug::%s]  at utg%.6d dir %d\n", __func__, (int)(vu>>1)+1, (int)(vu&1));}
        if (vu==(sink^1)){
            if (is_hinted_circular && (s.n==0)){  // special case: 1st pop of a hinted-circular subgraph, proceed to check child nodes
                ;
            }else{  // regular case: found end, end the walk
                stacku32_push(&s, vu);  // put back
                if (verbose) {fprintf(stderr, "[debug::%s]  reached sink\n", __func__);}
                break;
            }
        }
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        is_exhausted = 1;
        is_loop = 0;

        // try to greedily get the path with most coverage
        // so here we first sort the available target vertices
        //   before proceeding to topo checks
        //   (the sortting)
        stacku64_reset(&for_sorting);
        for (int i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            packed_value = ((uint64_t)ug->utg_coverage[av[i].v>>1])<<32 | (uint64_t)av[i].v;
            stacku64_push(&for_sorting, packed_value);
            if (verbose){
                fprintf(stderr, "[debug::%s]  (sort push) vertex %.6d , cov %.d\n", 
                    __func__, (int)(av[i].v>>1), (int)ug->utg_coverage[av[i].v>>1]);
            }
        }
        radix_sort_ovhamt64(for_sorting.a, for_sorting.a+for_sorting.n);
        stacku64_invert_buffer(&for_sorting);

        //   (topo checks)
        // for (int i=nv-1; i>=0; i--){ // TODO: less arbitrary
        for (int i=0; i<for_sorting.n; i++){  // greedily best
            wu = (uint32_t)for_sorting.a[i];
            // wu = av[i].v;
            if (stack32_is_in_stack(&visited, wu)) {
                if (wu==(sink^1)){
                    if (verbose) {fprintf(stderr, "[debug::%s]  reached sink(a)\n", __func__);}
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    is_reached_sink = 1;
                    break;
                }
                if (verbose) fprintf(stderr, "[debug::%s    (utg%.6d already visited)\n", __func__, (int)(wu>>1)+1);
                continue;
            }
            if (stack32_is_in_stack(&visited, wu^1)) {
                if (wu!=(sink^1)){
                    is_loop = 1;
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s]  reached sink(b)\n", __func__);}
                    is_reached_sink = 1;
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    break;
                }
            }
            if (verbose) fprintf(stderr, "[debug::%s]     not exhausted: %.6d dir %d\n", __func__, 
                                (int)(wu>>1)+1, (int)(wu&1));
            is_exhausted = 0;
            break;
        }

        if (is_reached_sink){break;}
        if (is_exhausted && is_loop){  // remove current node
            if (verbose) {fprintf(stderr, "[debug::%s]     found looping(a) at utg%.6d, remove %.6d\n", 
                                    __func__, (int)(wu>>1)+1, (int)(vu>>1)+1);}
                continue;
        }
        if (!is_exhausted){  // step
            if (stack32_is_in_stack(&visited, wu^1)) {  // check if it's a loop
                if (wu!=(sink^1)){
                    if (verbose) {fprintf(stderr, "[debug::%s]     found looping(b) at utg%.6d, remove %.6d\n", 
                                    __func__, (int)(wu>>1)+1, (int)(vu>>1)+1);}
                    continue;
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s]  reached sink(c)\n", __func__);}
                    is_reached_sink = 1;
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    break;
                }
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s]      push %.6d dir %d\n", __func__, (int)(wu>>1)+1, (int)(wu&1));}
                stacku32_push(&s, vu);  // put back
                stacku32_push(&s, wu);
                stacku32_push(&visited, wu);
            }
        }else{  // remove current node
            if (verbose) {fprintf(stderr, "[debug::%s]      removed %.6d\n", __func__, (int)(vu>>1)+1);}
            ;
        }
        san++;
        if (san>1000){
            fprintf(stderr, "[E::%s] too many steps\n", __func__);
            ret = 0;
            goto finish;
        }
    }
    if (s.n==0){
        if (verbose) {fprintf(stderr, "[E::%s] can't find a path\n", __func__);}
        ret = 0;
        goto finish;
    }else{
        if (verbose){
            if (verbose) {fprintf(stderr, "[debug::%s]  the path (%d nodes):\n", __func__, s.n);}
            for (int i=0; i<s.n; i++){
                fprintf(stderr, "        utg%.6d\n", (int)(s.a[i]>>1)+1);
            }
        }
        // hamt_ug_arc_flip_between(sg, ug, source, sink^1, 0, alt_label);  // mark vertices as alt
        stacku32_pop(&s, &vu);
        // recover the stored path
        while (stacku32_pop(&s, &wu)){
            // update black list
            if (wu!=source && wu!=sink){
                buf_blacklist[wu>>1] = 1;
            }
            
            hamt_ug_arc_del_selfcircle(sg, ug, wu, base_label);
            
            hamt_ug_arc_del(sg, ug, wu, vu, 0);
            // hamt_ug_utg_softdel(sg, ug, wu, base_label);
            // hamt_ug_utg_softdel(sg, ug, vu, base_label);
            nv = asg_arc_n(auxsg, wu);
            av = asg_arc_a(auxsg, wu);
            for (uint32_t i=0; i<nv; i++){
                if (av[i].v!=vu) {
                    color[av[i].v>>1] = 1;
                    hamt_ug_arc_del(sg, ug, wu, av[i].v, 1);
                    // hamt_ug_utg_softdel(sg, ug, av[i].v, alt_label);
                }else{
                    uu = vu^1;
                    nu = asg_arc_n(auxsg, uu);
                    au = asg_arc_a(auxsg, uu);
                    for (uint32_t j=0; j<nu; j++){
                        if (au[j].v!=(wu^1)){
                            hamt_ug_arc_del(sg, ug, uu, au[j].v, 1);
                            // hamt_ug_utg_softdel(sg, ug, av[i].v, alt_label);
                        }
                    }
                }
            }
            vu = wu;
        }

        ret = 1;
    }

finish:
    stacku32_destroy(&s);
    stacku32_destroy(&visited);
    stacku64_destroy(&for_sorting);
    free(color);
    return ret;
}
int hamt_ug_popTangles_v2(asg_t *sg, ma_ug_t *ug, uint32_t source, uint32_t sink, 
                       uint32_t *buf_blacklist,  // redundant but kept to be consistent with v1 interface
                       int base_label, int alt_label){
    // NOTE
    //     direction is source->...<-sink  // v1 interface consistency
    // RETURN
    //     number of arcs dropped
    int ret = hamt_ug_pop_subgraph(sg, ug, source, sink^1, base_label, alt_label);
    return ret;
}

int hamt_ug_resolveTangles(asg_t *sg, ma_ug_t *ug, 
                           int base_label, int alt_label){
    // FUNC
    //    Check and pop tangles. Heuristic.

    int nb_treated = 0, verbose = 0, max_label;
    int subgraph_size, non_tip_tigs;
    asg_t *auxsg = ug->g;
    uint32_t source, sink;  // source->...<-sink
    int max_length, max_coverage, has_long_tig;

    vecu32_t buf;
    vecu32_init(&buf);

    uint32_t *treated = (uint32_t*)calloc(auxsg->n_seq, sizeof(uint32_t));

    max_label = hamt_ug_util_BFS_markSubgraph(ug, base_label);
    if (verbose){fprintf(stderr, "[debug::%s] got %d subgraphs\n", __func__, max_label);}
    
    for (int idx_subgraph=0; idx_subgraph<max_label; idx_subgraph++){
        vecu32_reset(&buf);

        // if the subgraph has only a few tigs, or non of the tigs were long, skip it
        has_long_tig = 0;
        subgraph_size = non_tip_tigs = 0;
        for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
            // if (vu&1){continue;}  // only need to check  for one direction
            if (ug->u.a[vu>>1].len>500000) {has_long_tig = 1;}
            if (ug->u.a[vu>>1].subg_label==idx_subgraph){
                subgraph_size++;
                vecu32_push(&buf, vu);
            }
            if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 &&
                hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=0){
                    non_tip_tigs++;
            }
        }
        if (subgraph_size<2 || non_tip_tigs<2 || !has_long_tig){
            continue;
        }

        // check pairwise
        for (int idx_vu=0; idx_vu<buf.n; idx_vu++){
            for (int idx_wu=0; idx_wu<buf.n; idx_wu++){
                max_length = max_coverage = 0;
                source = buf.a[idx_vu];
                sink = buf.a[idx_wu]^1;
                if (treated[source>>1] || treated[sink>>1]){
                    if (verbose) {fprintf(stderr, "[debug::%s] skip becuase either handle was a part of some tangles treated\n", __func__);}
                    continue;
                }

                if (verbose){
                    fprintf(stderr, "[debug::%s] check source %.6d sink %.6d \n", __func__, 
                                            (int)(source>>1)+1, (int)(sink>>1)+1);
                }

                // // debug
                // if (((source>>1)+1)!=1374) {continue;}
                // if (((sink>>1)+1)!=2694) {continue;}

                // (note: must check label, since hamt_ug_arc_flip_between is used and it only flip labels for simplicity)
                if (base_label>=0 && auxsg->seq_vis[source>>1]!=base_label) {continue;}
                if (base_label>=0 && auxsg->seq_vis[sink>>1]!=base_label) {continue;}

                if (source==sink) {  // for hinted circular, it's source==(sink^1)
                    if (verbose) fprintf(stderr, "[deubg::%s] skip becuase sink==source\n", __func__);
                    continue;
                } 
                if (hamt_ug_check_localSourceSinkPair(ug, source, sink, &max_length, &max_coverage, base_label, asm_opt.gc_tangle_max_tig)){
                    if ( ((ug->u.a[source>>1].len<50000) && (ug->u.a[sink>>1].len<50000)) ||   // note/TODO/bug: this is an arbitrary waterproof, trying to prevent weird cases (since we searched both directions during traversal to tolerate inversions/loops)
                         ( (max_length>500000)) ||
                         (max_coverage > (ug->utg_coverage[source>>1]*2)) ||
                         (max_coverage > (ug->utg_coverage[sink>>1]  *2))
                         ){
                        if (verbose){
                            fprintf(stderr, "[debug::%s] failed length or cov. Max l %d, max cov %d (handles' are %d and %d)\n", 
                                                __func__, 
                                                max_length, max_coverage, (int)(ug->utg_coverage[source>>1]), (int)ug->utg_coverage[sink>>1]);
                        }
                    }else{
                        if (verbose){fprintf(stderr, "[debug::%s] try to pop tangle\n", __func__);}
                        if (hamt_ug_popTangles(sg, ug, source, sink, treated, base_label, alt_label)){
                            nb_treated++;
                        }else{
                            if ((source>>1)!=(sink>>1)){
                                if (verbose) {fprintf(stderr, "[W::%s] did not pop: utg%.6d - utg%.6d\n", __func__, (int)(source>>1)+1, (int)(sink>>1)+1);}
                            }
                        }
                    }
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s] source-sink pair check negative\n", __func__);}
                }
            }
        }

    }

    if (nb_treated){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        hamt_ug_cleanup_arc_by_labels(sg, ug);
    }

    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n", __func__, nb_treated);
    }
    vecu32_destroy(&buf);
    free(treated);
    return nb_treated;
}

int hamt_ug_resolve_fake_haplotype_bifurcation(asg_t *sg, ma_ug_t *ug, int base_label,
                                               ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources){
    // NOTE
    //    Naming isn't informative; a real case is the mock dataset, cutibacterium_acnes subgraph (tested in r10 and before),
    //      where the short unitigs in the center aren't haplotips to each other. The reads aligned, however, persumably due to
    //      large overhang, they were ignored and not logged in R_INF.paf and R_INF.reverse_paf (good haplotigs
    //      should have some overlaps logged in reverse_paf). Not entirely sure why this ends up with the asg&ug topo,
    //      but by examining paf, reverse_paf and OVEC_INF (which logs bad but calculated overlaps during the final ovec),
    //      at least for this particular case, the correct walk can be recovered 
    //    (heuristic: "correct" pathway has better connection before trans reduct). 
    //    (A self note for the above case: walk is readname 909159->1028753->5483449->2242924)
    // FUNC
    //    see code & comments, lump of experimental heuristics.
    // TODO
    //    currently only solving the seen case described above. relax the topo check? (e.g. instead of two ends of uu, let it be two contigs) 
    /* 
              wu1---uu_end1
             /            | \...
    -------vu             |
            \             |
             wu2----uu_end2-...
        wu1 and wu2 are short and doesn't appear to be haplotigs; vu and uu have direct overlap; uu may have other connections.
    */
    int verbose = 0;
    int nb_treated = 0;

    asg_t *auxsg = ug->g;
    uint32_t wu[2], uu[2], nv, u[2], v;
    asg_arc_t *av;
    ma_hit_t *h, *h_rev;
    int stat[2];
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=2){
            continue;
        }
        hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label);
        if (ug->u.a[wu[0]>>1].len>30000){continue;}
        if (ug->u.a[wu[1]>>1].len>30000){continue;}
        if (hamt_asgarc_util_countSuc(auxsg, wu[0], 0, 0, base_label)!=1 ||
            hamt_asgarc_util_countSuc(auxsg, wu[0], 0, 0, base_label)!=1){
                continue;
            }

        if (verbose) {fprintf(stderr, "[debug::%s] at utg%.6d, targets utg%.6d utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);}
        // check and get uu
        hamt_asgarc_util_get_the_one_target(auxsg, wu[0], &uu[0], 0, 0, base_label);
        hamt_asgarc_util_get_the_one_target(auxsg, wu[1], &uu[1], 0, 0, base_label);
        if (uu[0]!=(uu[1]^1)){
            if (verbose) {fprintf(stderr, "[debug::%s]     uu didn't pass\n", __func__);}
            continue;
        }  // note: it's ok for uu to have backward branching, so not checking that.

        // check overlap status
        if (hamt_check_diploid(ug, wu[0], wu[1], 0.3, reverse_sources)<=0 && 
            hamt_check_suspicious_diploid(sg, ug, wu[0], wu[1], 0.3)>0){
            // ^i.e. wu1 and wu2 are not reliable haplotigs, but they also have a decent amount of not-that-good overlaps
            // check: vu vs uu
            if (uu[0]&1){
                u[0] = ug->u.a[uu[0]>>1].end;
                u[1] = ug->u.a[uu[0]>>1].start;
            }else{
                u[0] = ug->u.a[uu[0]>>1].start;
                u[1] = ug->u.a[uu[0]>>1].end;
            }
            if (vu&1){
                v = ug->u.a[vu>>1].start^1;
            }else{
                v = ug->u.a[vu>>1].end^1;
            }
            stat[0] = does_ovlp_ever_exist(sources, v, u[0], &h, &h_rev);
            stat[1] = does_ovlp_ever_exist(sources, v, u[1], &h, &h_rev);
            if (stat[0] && stat[1]){  // if two looked both good, do nothing
                if (verbose) {fprintf(stderr, "[debug::%s]     vu overlaps with both ends of uu\n", __func__);}
                continue;
            } 
            if (stat[0]){
                hamt_ug_arc_del(sg, ug, vu, wu[1], 1);
                hamt_ug_arc_del(sg, ug, wu[1], uu[1],1);
                hamt_ug_utg_softdel(sg, ug, wu[1], 1);
                if (verbose) {fprintf(stderr, "[debug::%s]     dropped utg%.6d\n", __func__, (int)(wu[1]>>1)+1);}
                nb_treated++;
            }else if (stat[1]){
                hamt_ug_arc_del(sg, ug, vu, wu[0], 1);
                hamt_ug_arc_del(sg, ug, wu[0], uu[0],1);
                hamt_ug_utg_softdel(sg, ug, wu[0], 1);
                nb_treated++;
                if (verbose) {fprintf(stderr, "[debug::%s]     dropped utg%.6d\n", __func__, (int)(wu[0]>>1)+1);}
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s]     vu vs uu failed\n", __func__);}
            }
        }else{
            if (verbose) {fprintf(stderr, "[debug::%s]     dip check didn't pass\n", __func__);}
        }
    }

    if (nb_treated){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n\n", __func__, nb_treated);
    }
    return nb_treated;
}

int hamt_ug_resolve_fake_haplotype_bifurcation_aggressive(asg_t *sg, ma_ug_t *ug, int base_label,
                                               ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources){
    // NOTE
    //     hamt_ug_resolve_fake_haplotype_bifurcation but doesn't require the uu (i.e. wu1 and wu2 can target different unitigs)
    //     Will check vu vs wu1's target and wu2's target more rigorously than the original function,
    //       since the very start/last reads might not form a good overlap (e.g. "containment" with long overhang)
    /* 
          ...
             \
              wu1---uu1...
             /   \---...
    -------vu
            \   /---...
             wu2----uu2...
            /
         ...
        wu1 and wu2 are short and doesn't appear to be haplotigs; 
        vu and [uu1 or uu2] to have direct overlap; 
        uu1 and uu2 may have other connections.
        wu1 and wu2 may have other backward connections, 
        vu does NOT have other backward connections, just to be safe.
    // FUNC / NOTE
        Essentially, this fucntion means that since we expect wu to be very short,
          if a walk is a "good" or "normal" one, vu is very likely to have some acceptable overlap(s)
          directly with uu. If wu1 and wu2 doesn't appear to be a pair of haplotig, then we want to check
          this. And then if one of them fulfills this expectation but not the other, we may want to ditch the later.
        Generalized from hamt_ug_resolve_fake_haplotype_bifurcation, this is a more dangerous version.
          If more data suggests against it, this should be taken off. (potential bug / TODO)
    */
    int verbose = 0;
    int nb_treated = 0;

    asg_t *auxsg = ug->g;
    uint32_t wu[2], uu[2], nv, u1, u2, v;
    asg_arc_t *av;
    ma_hit_t *h, *h_rev;
    int stat[2];
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=2){
            continue;
        }
        hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label);
        if (ug->u.a[wu[0]>>1].len>30000){continue;}
        if (ug->u.a[wu[1]>>1].len>30000){continue;}
        if (hamt_asgarc_util_countSuc(auxsg, wu[0], 0, 0, base_label)!=1 ||
            hamt_asgarc_util_countSuc(auxsg, wu[0], 0, 0, base_label)!=1){
                continue;
            }

        if (verbose) {fprintf(stderr, "[debug::%s] at utg%.6d, targets utg%.6d utg%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);}
        // check and get uu
        hamt_asgarc_util_get_the_one_target(auxsg, wu[0], &uu[0], 0, 0, base_label);
        hamt_asgarc_util_get_the_one_target(auxsg, wu[1], &uu[1], 0, 0, base_label);
        if ((uu[0]>>1)==(uu[1]>>1)){continue;}  // handle these with other routines

        // check overlap status
        if (hamt_check_diploid(ug, wu[0], wu[1], 0.3, reverse_sources)<=0 && 
            hamt_check_suspicious_diploid(sg, ug, wu[0], wu[1], 0.3)>0){
            // ^i.e. wu1 and wu2 are not reliable haplotigs, but they also have a decent amount of not-that-good overlaps
            // check: vu vs uu
            if (uu[0]&1){
                u1 = ug->u.a[uu[0]>>1].end;
            }else{
                u1 = ug->u.a[uu[0]>>1].start;
            }
            if (uu[1]&1){
                u2 = ug->u.a[uu[1]>>1].end;
            }else{
                u2 = ug->u.a[uu[1]>>1].start;
            }
            if (vu&1){
                v = ug->u.a[vu>>1].start^1;
            }else{
                v = ug->u.a[vu>>1].end^1;
            }
            stat[0] = does_ovlp_ever_exist(sources, v, u1, &h, &h_rev);
            stat[1] = does_ovlp_ever_exist(sources, v, u2, &h, &h_rev);
            if (stat[0] && stat[1]){  // if two looked both good, do nothing
                if (verbose) {fprintf(stderr, "[debug::%s]     vu overlaps with both ends of uu\n", __func__);}
                continue;
            } 
            if (stat[0]){
                hamt_ug_arc_del(sg, ug, vu, wu[1], 1);
                hamt_ug_arc_del(sg, ug, wu[1], uu[1],1);
                hamt_ug_utg_softdel(sg, ug, wu[1], 1);
                if (verbose) {fprintf(stderr, "[debug::%s]     dropped utg%.6d\n", __func__, (int)(wu[1]>>1)+1);}
                nb_treated++;
            }else if (stat[1]){
                hamt_ug_arc_del(sg, ug, vu, wu[0], 1);
                hamt_ug_arc_del(sg, ug, wu[0], uu[0],1);
                hamt_ug_utg_softdel(sg, ug, wu[0], 1);
                nb_treated++;
                if (verbose) {fprintf(stderr, "[debug::%s]     dropped utg%.6d\n", __func__, (int)(wu[0]>>1)+1);}
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s]     vu vs uu failed\n", __func__);}
            }
        }else{
            if (verbose) {fprintf(stderr, "[debug::%s]     dip check didn't pass\n", __func__);}
        }
    }

    if (nb_treated){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] treated %d spots\n\n", __func__, nb_treated);
    }
    return nb_treated;
}

void hamt_debug_get_diploid_info_about_all_branchings(ma_ug_t *ug, ma_hit_t_alloc *reverse_sources){
    float ratio1, ratio2;
    asg_t *auxsg = ug->g;
    vecu32_t v;
    vecu32_init(&v);
    uint32_t nv;
    asg_arc_t *av;
    int alarm;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        vecu32_reset(&v);
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)<2){continue;}
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);    
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, 0)>0){
                vecu32_push(&v, av[i].v);
            }
        }
        if (v.n<=1){continue;}
        // check diploid, pairwise
        fprintf(stderr, "[debug::%s] at utg%.6d\n",__func__, (int)(vu>>1)+1);
        alarm = 1;
        for (int i=0; i<v.n-1; i++){
            for (int j=i+1; j<v.n; j++){
                // het overlap
                hamt_check_diploid_report(ug, v.a[i], v.a[j], reverse_sources, &ratio1, &ratio2);
                fprintf(stderr, "[debug::%s]    utg%.6d vs utg%.6d, ratios: %.2f and %.2f\n", 
                                            __func__,
                                            (int)(v.a[i]>>1)+1, (v.a[j]>>1)+1,
                                            ratio1, ratio2);
                if (ratio1>0 || ratio2>0) {alarm = 0;}
                
                // suspicious overlap
                hamt_check_suspicious_diploid_report(NULL, ug, v.a[i], v.a[j], &ratio1, &ratio2);
                fprintf(stderr, "[debug::%s]    utg%.6d vs utg%.6d, sus ratios: %.2f and %.2f\n", 
                                            __func__,
                                            (int)(v.a[i]>>1)+1, (v.a[j]>>1)+1,
                                            ratio1, ratio2);
            }
        }
        if (alarm){
            fprintf(stderr, "[debug::%s] !!! none of the pairs check check out !!!\n",__func__);
        }
    }
    vecu32_destroy(&v);
}

int hamt_ug_cleanup_almost_circular(asg_t *sg, ma_ug_t *ug, int base_label){
    // FUNC
    //     If a long contig forms a circle, drop tips if the tip is relatively short (abs length can be huge)
    //     (This is seen in moderate coverage cases (~x30), the tip might be a broken haplotype?)
    int verbose = 0;
    int ret = 0, checksout=0;
    int contig_length_threshold=3000000;
    float ratio = 0.3;

    asg_t *auxsg = ug->g;
    uint32_t vu, nv;
    asg_arc_t *av;

    // drop tips only
    for (vu=0; vu<auxsg->n_seq; vu++){
        checksout = 0;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            if (av[i].v==vu){
                checksout = 1;
                break;
            }
        }
        if (!checksout) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==1) continue;
        if (ug->u.a[vu>>1].len<contig_length_threshold) continue;

        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            if (av[i].v==vu) continue;
            if (hamt_asgarc_util_isTip(auxsg, av[i].v, 0, 0, base_label)){
                hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
                ret+=1;
                if (verbose){
                    fprintf(stderr, "[debug::%s] base utg%.6d, dropped tip utg%.6d\n", __func__, 
                                        (int)(vu>>1)+1, (int)(av[i].v>>1)+1);
                }
            }
        }
    }

    if (ret){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return ret;
}

int hamt_ug_rescue_bifurTip(asg_t *sg, ma_ug_t *ug, int base_label,
                           ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, const ma_sub_t* coverage_cut){
    // NOTE/BUG/TODO
    //     The checking will only check if there is any homo overlaps between
    //      the two targets; however, recover requires the end reads to 
    //      overlap with each other.
    //     This might be problemetic in weird situations.
    int verbose = 0;
    asg_t *auxsg = ug->g;
    int cnt = 0;

    uint32_t wu[2];
    uint32_t startread, endread;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len>50000) continue;  // handle utg shall be short

        if (hamt_ug_check_bifurTip(ug, vu, base_label, sources, reverse_sources)>0){
            hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label);

            // get reads at the ends of wu0 and wu1
            if (wu[0]^1){
                startread = ug->u.a[wu[0]>>1].end^1;
            }else{
                startread = ug->u.a[wu[0]>>1].start^1;
            }
            if (wu[1]^1){
                endread = ug->u.a[wu[1]>>1].end;
            }else{
                endread = ug->u.a[wu[1]>>1].start;
            }

            
            if (hamt_ug_recover_ovlp_if_existed(sg, ug, startread, endread, sources, coverage_cut, 0, 1)>0){
                hamt_ug_arc_del(sg, ug, vu, wu[0], 1);
                hamt_ug_arc_del(sg, ug, vu, wu[1], 1);
                cnt+=1;
                if (verbose){
                    fprintf(stderr, "[debug::%s] treated, handle %.6d side1 %.6d side2 %.6d\n", __func__, 
                                        (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);
                }
            }else{
                if (verbose){
                    fprintf(stderr, "[debug::%s] did NOT treat bcs can't recover arc, handle %.6d side1 %.6d side2 %.6d\n", __func__, 
                                        (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);
                }
            }

            

        }
    }


    if (cnt){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return cnt;
}

// DO NOT USE (maybe)
int hamt_asg_arc_del_intersample_branching(asg_t *sg,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex){
    // FUNC
    //     To be called after transitive reduction and basic asg cleaning routines.
    //     For example, if we see a branching at segment A (sample X) whose child segments are 
    //       B (sample X) and C (sample Y), we might want to favor B and discard A-C connection.
    // NOTE
    //     Will generate a temp ug.
    // RET
    //     Number of arcs dropped.

    assert(asm_opt.mode_coasm);
    if (!R_INF.coasm_sampleID){
        fprintf(stderr, "[E::%s] co-assembly, sample info hasn't been loaded?\n", __func__);
        exit(1);
    }

    int ret = 0;
    int verbose = 0;
    double startTime = Get_T();
    ma_ug_t *ug = ma_ug_gen(sg);
    asg_t *auxsg = ug->g;

    uint32_t nv, wu; 
    asg_arc_t *av;
    uint32_t handle, target[2], target_v[2];

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countNoneTipSuc(auxsg, vu, -1)!=2){continue;}
        if (ug->utg_coverage[vu>>1]>10) {continue;}
        handle = ug->u.a[vu>>1].end;  // seqID, with direction bit
        
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        int san_i = 0;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) {continue;}
            wu = av[i].v;
            if (hamt_asgarc_util_isTip(auxsg, wu, 0, 0, -1)) {continue;}
            target[san_i] = ug->u.a[wu>>1].start;
            target_v[san_i] = wu;
            san_i++;
        }
        if (san_i!=2){
            fprintf(stderr, "[E::%s] should not happen, san_i is %d\n", __func__, san_i);
            exit(1);
        }

        // check sampleIDs of the segments
        if (R_INF.coasm_sampleID[target[0]>>1]==R_INF.coasm_sampleID[target[1]>>1]){
            continue;
        }
        if (R_INF.coasm_sampleID[handle>>1]!=R_INF.coasm_sampleID[target[0]>>1] &&
            R_INF.coasm_sampleID[handle>>1]!=R_INF.coasm_sampleID[target[1]>>1]){
            continue;
        }
        if (R_INF.coasm_sampleID[handle>>1]==R_INF.coasm_sampleID[target[0]>>1]){
            hamt_ug_arc_del(sg, ug, vu, target_v[1], 1);
        }else{
            hamt_ug_arc_del(sg, ug, vu, target_v[0], 1);
        }
        ret+=1;
    }
    if (ret){
        asg_cleanup(sg);
        // asg_cleanup(auxsg);  // no need
    }

    hamt_ug_destroy(ug);
    // (clean up sg)
    if (verbose){
        fprintf(stderr, "[debug::%s] treated %d spots, used %.2f s.\n", __func__, ret, Get_T()-startTime);
    }
    return ret;
}

int hamt_ug_cut_very_short_multi_tip(asg_t *sg, ma_ug_t *ug, int nb_threshold){
    // FUNC
    //     drop multi-tip if it's less than nb_threshold reads
    // RET
    //     number of spots treated
    int ret = 0;
    asg_t *auxsg = ug->g;
    uint32_t nv;
    asg_arc_t *av;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, -1)!=0) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, -1)<2) continue;
        if (ug->u.a[vu>>1].n>nb_threshold) continue;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
        }
        ug->u.a[vu>>1].c = 1;
        ug->g->seq_vis[vu>>1] = 1;
        ret++;
    }
    if (ret){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return ret;
}



int hamt_ug_drop_shorter_ovlp(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources){
    // NOTE
    //     Assume coverage has been collected.
    // FUNC
    //     Drop the shorter ovlp at a birfucation if the length diff between the two 
    //       is extreme. Additionally, will consider:
    //       - coverage diff
    //       - topo: query has at least 2 non-tip targets (won't check further though)
    // RETURN
    //     number of ovlps dropped.
    int verbose = 0;
    int ret = 0;
    asg_t *auxsg = ug->g;
    int buf_l;
    uint32_t nv, wu, w, v;  // w is the start of wu, v is the end^1 of vu
    asg_arc_t *av;
    
    ma_hit_t *h;
    ma_hit_t *h_tmp;
    uint32_t h_vu, h_wu;
    int worst_cov, cov, max_cov;
    int worst_ovlp, ovlpl, max_ovlp;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (auxsg->seq_vis[vu>>1]!=0) {continue;}
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)<2){continue;}

        // reset
        buf_l = 0;
        worst_cov = 1000;  // a large number
        worst_ovlp = 100000;  // a large number
        max_cov = 0;
        max_ovlp = 0;

        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) continue;
            if (auxsg->seq_vis[av[i].v>>1]!=0) continue;
            if (hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, 0)==0) continue;  // target is a tip
            wu = av[i].v;
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
            
            h_tmp = get_specific_overlap_handle(sources, v>>1, w>>1);
            if (!h_tmp) {  // can't get the handle for some reason
                h_tmp = get_specific_overlap_handle(reverse_sources, v>>1, w>>1);
                if (!h_tmp){
                    fprintf(stderr, "[E::%s] can't get handle for v %.*s (%d) w %.*s (%d) (no dir bit)\n",
                                        __func__, 
                                        (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)(v>>1), 
                                        (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1), (int)(w>>1));
                    continue;
                }
            }
            if (!h_tmp->el) continue;  // skip inexact overlaps

            cov = ug->utg_coverage[wu>>1];
            ovlpl = h_tmp->te-h_tmp->ts;
            max_ovlp = ovlpl>max_ovlp? ovlpl : max_ovlp;
            max_cov = cov>max_cov? cov : max_cov;
            
            if ((float)ovlpl/max_ovlp > 0.7) continue;  // diff not significant
            if ((float)cov/max_cov>0.7) continue;  // diff not significant
            if (ovlpl<worst_ovlp && cov<(float)worst_cov*1.5){
                h = h_tmp;
                worst_ovlp = ovlpl;
                worst_cov = cov;
                h_vu = vu;
                h_wu = wu;
                buf_l++;
            }
        }
        if (buf_l==0) continue;  // nothing on the radar

        // if we reach here, drop stuff
        hamt_ug_arc_del(sg, ug, h_vu, h_wu, 1);
        ret++;
        if (verbose){
            fprintf(stderr, "[debug::%s] cut between %.*s and %.*s \n",
                                __func__, 
                                (int)Get_NAME_LENGTH(R_INF, h_vu>>1), Get_NAME(R_INF, h_vu>>1),
                                (int)Get_NAME_LENGTH(R_INF, h_wu>>1), Get_NAME(R_INF, h_wu>>1));
            fprintf(stderr, "[debug::%s]    max_cov %d worst_cov %d max_ov %d worst_ov %d\n", 
                                __func__,
                                max_cov, worst_cov, max_ovlp, worst_ovlp);
        }
    }
    fprintf(stderr, "[M::%s] cut %d\n", __func__, ret);
    return ret;
}


// int hamt_ug_drop_worse_ovlp_at_bifur(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc *sources){
//     // NOTE
//     //      Was written during r43->r44, dropped because worse empirical performance and 
//     //       potentail overfitting.            
//     // FUNC
//     //     For example, given a->b (inexact) and a->c (exact), 
//     //      where a has 18x coverage, b 34x, c 16x,
//     //      we might want to drop a->b.
//     // RETURN
//     //     Number of arcs dropped.
//     int ret = 0;
//     asg_t *auxsg = ug->g;
//     ma_hit_t *h;

//     uint32_t vu, nv, wu, v, w;
//     asg_arc_t *av;
//     int sancheck;
//     int is_exact;
//     uint32_t targets[5];
//     int covs[5];
//     int ovlp_lengths[5];
//     int exact_status[5];

//     for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
//         if (hamt_asgarc_util_countNoneDanglingTipSuc(auxsg, vu, 0)!=2){continue;}
//         nv = asg_arc_n(auxsg, vu);
//         av = asg_arc_a(auxsg, vu);
        
//         sancheck = 0;
//         is_exact = 0;
//         for (int i=0; i<nv; i++){
//             if (av[i].del) continue;
//             wu = av[i].v;
//             if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, 0)==0 && 
//                 hamt_asgarc_util_countPre(auxsg, wu, 0, 0, 0)==1) continue;  // is a tip with only the vu as its neighbor

//             if ((vu&1)==0){
//                 v = ug->u.a[vu>>1].end^1;
//             }else{
//                 v = ug->u.a[vu>>1].start^1;
//             }
//             if ((wu&1)==0){
//                 w = ug->u.a[wu>>1].start;
//             }else{
//                 w = ug->u.a[wu>>1].end;
//             }

//             h = get_specific_overlap_handle(sources, v>>1, w>>1);
//             if (!h) continue;  // ERROR
//             if (h->el) is_exact++;
            
//             targets[sancheck] = wu;
//             covs[sancheck] = HAMT_DIFF(ug->utg_coverage[wu>>1], ug->utg_coverage[vu>>1]);
//             ovlp_lengths[sancheck] = (int)(h->te - h->ts);
//             exact_status[sancheck] = (int)h->el;
//             sancheck++;
//             if (sancheck>2) break;
//         }
//         // assert(sancheck==2);
//         if (sancheck!=2) continue; // ERROR
        
//         if (is_exact!=1) continue;  // either both are inexact, or both are exact. Don't try to treat these.
//         if (HAMT_MIN(covs[0], covs[1])>10) continue;  // one of the target contig shall have similar coverage to the query contig
//         if ((float)HAMT_MAX(covs[0], covs[1])/HAMT_MIN(covs[0], covs[1])<5) continue;  // ... and the other target does not.

//         // if we reach here, drop the worse one
//         if (exact_status[0]==0) {
//             wu = targets[1];
//         }else{
//             wu = targets[0];
//         }
//         hamt_ug_arc_del(sg, ug, vu, wu, 1);
//         ret++;
        
//     }


//     fprintf(stderr, "[M::%s] dropped %d\n", __func__, ret);
//     return ret;
// }







int hamt_ug_3mer_cut_with_coverage(asg_t *sg, ma_ug_t *ug, int threshold_l){
    // FUNC
    //     This is meant to be part of the final pruning.
    //     Collect all trinucleotide profile of long contigs 
    //      (also assume that unitig coverage has been made available),
    //      for each bifurcation where all 3 unitigs involved are long, 
    //      pick one path if trinucleotide profile and coverage hints agree.
    // TODO
    //     Waterproof the math; in meta assembly we don't have super long contigs right now though.
    // RETURN
    //     Bifurcations treated.
    int verbose = 0;
    int ret = 0;
    asg_t *auxsg = ug->g;
    uint32_t *map = (uint32_t*)calloc(auxsg->n_seq, sizeof(uint32_t));  // unitig ID to indices used by the linear buffer
    float **profile;
    float **profile_rev;
    uint32_t idx = 0;
    ma_utg_t *p;
    char *seq;

    // count long unitigs
    for (uint32_t vu=0; vu<auxsg->n_seq; vu++){
        map[vu] = idx;
        if (ug->u.a[vu].len>=__FLT_MAX__){
            fprintf(stderr, "[W::%s] contig too long, skipping it\n", __func__);
        }
        if (ug->u.a[vu].len>threshold_l){
            idx++;
        }
    }
    if (verbose){fprintf(stderr, "[debug::%s] total %d long tigs\n", __func__, (int)idx);}
    profile = (float**)calloc(idx, sizeof(float*));
    profile_rev = (float**)calloc(idx, sizeof(float*));
    for (int i=0; i<idx; i++){
        profile[i] = (float*)calloc(64, sizeof(float));
        profile_rev[i] = (float*)calloc(64, sizeof(float));
    }

    // collect 3mer profiles, assuming no N base
    idx = 0;
    uint32_t tmp;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) continue;
        if (ug->u.a[vu>>1].len<=threshold_l) continue;
        p = &ug->u.a[vu>>1];
        seq = p->s;
        assert(seq);
        for (int i=0; i<ug->u.a[vu>>1].len-3+1; i++){
            tmp = seq_nt4_table[(uint8_t)seq[i]]*16 + seq_nt4_table[(uint8_t)seq[i+1]]*4 + seq_nt4_table[(uint8_t)seq[i+2]];
            profile[idx][tmp]++;
            tmp = (3-seq_nt4_table[(uint8_t)seq[i+2]])*16 + (3-seq_nt4_table[(uint8_t)seq[i+1]])*4 + (3-seq_nt4_table[(uint8_t)seq[i]]);
            profile_rev[idx][tmp]++;
        }
        for (int i=0; i<64; i++){
            profile[idx][i] = profile[idx][i]/(ug->u.a[vu>>1].len-2);
            profile_rev[idx][i] = profile_rev[idx][i]/(ug->u.a[vu>>1].len-2);
        }
        idx++;
    }

    // check bifurcations
    uint32_t wu[2];
    float cos1, cos2;
    float *h_vu, *h_wu1, *h_wu2;
    int covdiff1, covdiff2;
    int which_cos, which_cov;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len<=threshold_l || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)!=2) continue;
        hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, 0);
        if (ug->u.a[wu[0]>>1].len<=threshold_l || ug->u.a[wu[1]>>1].len<=threshold_l) continue;  // require both targets to be long
        // if (hamt_asgarc_util_countSuc(auxsg, wu[0], 0, 0, 0)==0 || hamt_asgarc_util_countSuc(auxsg, wu[1], 0, 0, 0)==0) continue;  // targets are not tips
        
        // get profile hints
        if (verbose){fprintf(stderr, "[debug::%s] > %.6d\n", __func__, (int)(vu>>1)+1);}
        if (verbose){fprintf(stderr, "[debug::%s]     get profile\n", __func__);}
        if(vu&1){h_vu = profile_rev[map[vu>>1]];}
        else{h_vu = profile[map[vu>>1]];}
        if (wu[0]&1){h_wu1 = profile_rev[map[wu[0]>>1]];}
        else{h_wu1 = profile[map[wu[0]>>1]];}
        if (wu[1]&1){h_wu2 = profile_rev[map[wu[1]>>1]];}
        else{h_wu2 = profile[map[wu[1]>>1]];}
        cos1 = 1-cosine_similarity(h_vu, h_wu1, 64);
        cos2 = 1-cosine_similarity(h_vu, h_wu2, 64);
        if (cos1>0.01 && cos2>0.01) {  // large cosine distances - somehow both ways look like bad choices
            if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0){  // handle is a tip, cut both
                hamt_ug_arc_del(sg, ug, vu, wu[0], 1);
                hamt_ug_arc_del(sg, ug, vu, wu[1], 1);
                if (VERBOSE || verbose){
                    fprintf(stderr, "[M::%s] ctg%.6d -> ctg%.6d/ctg%.6d, drop handle\n", __func__,
                                    (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);
                }
            }else{
                if (verbose){
                    fprintf(stderr, "[M::%s] ctg%.6d -> ctg%.6d/ctg%.6d, cosine does not support\n", __func__,
                                    (int)(vu>>1)+1, (int)(wu[0]>>1)+1, (int)(wu[1]>>1)+1);
                }
            }
            continue;
        }else{
            if (cos1-cos2>0.005){
                which_cos = 1;  // prefer wu[1]
            }else if (cos2-cos1>0.005){
                which_cos = 0;  // prefer wu[0]
            }else{
                if (verbose){fprintf(stderr, "[debug::%s]     profile diff not significant\n", __func__);}
                continue;  // profile diff not significant
            }
        }
        
        // get coverage hints
        if (verbose){fprintf(stderr, "[debug::%s]     get coverage\n", __func__);}
        covdiff1 = abs(ug->utg_coverage[wu[0]>>1] - ug->utg_coverage[vu>>1]);
        covdiff2 = abs(ug->utg_coverage[wu[1]>>1] - ug->utg_coverage[vu>>1]);
        if (abs(covdiff1-covdiff2)<5) {
            if (verbose){fprintf(stderr, "[debug::%s]     coverage diff not significant\n", __func__);}
            continue;
        }
        if (covdiff1<covdiff2) {which_cov = 0;}
        else {which_cov = 1;}

        // do they agree?
        if (which_cos==which_cov){
            if (VERBOSE || verbose){
                fprintf(stderr, "[M::%s] ctg%.6d -> ctg%.6d, discard path to ctg%.6d\n", __func__,
                                 (int)(vu>>1)+1, (int)(wu[which_cos]>>1)+1, (int)(wu[!which_cos]>>1)+1);
            }
            hamt_ug_arc_del(sg, ug, vu, wu[!which_cos], 1);
            ret++;
        }else{
            if (verbose){fprintf(stderr, "[debug::%s]     decision did not agree, pass\n", __func__);}
        }
    }

    free(map);
    for (int i=0; i<idx; i++){free(profile[i]); free(profile_rev[i]);}
    free(profile);
    free(profile_rev);
    // if (ret){
    //     asg_cleanup(sg);
    //     asg_cleanup(auxsg);
    // }
    if (verbose){fprintf(stderr, "[debug::%s] total %d\n", __func__, ret);}
    return ret;
}


int hamt_ug_finalprune(asg_t *sg, ma_ug_t *ug){
    // FUNC
    //     A set of aggressive pruning after p_ctg. 
    //     Will drop very long tips (>100kb) etc. Doing this because
    //      currently binning does not use or cannot efficiently use
    //      assembly graph's topology information, and in additional,
    //      might want to recruit contigs that are not *that* close
    //      to each other on the assembly graph.
    //     Therefore it might be of
    //      interest to just let hifiasm-meta get rid of dangling stuff
    //      and report longer contigs.
    //     Better binning or alike would be more preferable in the future.
    // TODO
    //     To duplicate or not to duplicate.
    // RETURN
    //     Number of treatments.
    int verbose = 0;
    int ret = 0;

    asg_t *auxsg = ug->g;

    uint32_t nv, vu, wu;
    asg_arc_t *av;

    // "binning" for long tigs 
    ret += hamt_ug_3mer_cut_with_coverage(sg, ug, 200000);
    fprintf(stderr, "[M::%s] 3mer cut finished, total treatments so far: %d\n", __func__, ret);

    // topo-dependent pruning of tips, regardless of the tig's length
    // TODO: should consider target lengths
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)!=0 || hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)!=1) continue;
        if (hamt_asgarc_util_get_the_one_target(auxsg, vu, &wu, 0, 0, 0)<0){
            fprintf(stderr, "[E::%s] can't get target for ctg%.6d when should\n", __func__, (int)(vu>>1)+1);
            continue;
        }

        if (verbose){fprintf(stderr, "[debug::%s] > %.6d\n", __func__, (int)(vu>>1)+1);}
        if (hamt_asgarc_util_countSuc(auxsg, wu, 0, 0, 0)==0 || hamt_asgarc_util_countNoneTipSuc(auxsg, wu^1, 0)==0) continue;  // target must has other links in both direction
        if (verbose){fprintf(stderr, "[debug::%s]     checkpoint\n", __func__);}

        hamt_ug_arc_del(sg, ug, vu, wu, 1);
        if (VERBOSE){
            fprintf(stderr, "[M::%s] long tig drop: ctg%.6d\n", __func__, (int)(vu>>1)+1);
        }
        ret++;
    }

    // topo-dependent pruning of small subgraphs 

    // clean
    if (ret){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    fprintf(stderr, "[M::%s] all finished, total treatments: %d\n", __func__, ret);
    return ret;
}


int hamt_ug_popLooseTangles(asg_t *sg, ma_ug_t *ug, int threshold_min_handle_length){
    // FUNC
    //      An aggressive version of tangle popping.
    //     `hamt_ug_popTangles` requires the tangle to be 
    //       enclosed by a source and a sink.
    //      This function relaxes this criteria and tries to
    //       form (partially phased) unitigs unless the path
    //       absolutely has to part ways.
    //      One example is the assembly of 4 E.coli strains (pileup zymo).
    //      In the force-directed graph layout by Bandage, the graph has
    //       visual "bubbles". The complex local structures on each of these
    //       "bubble edges" are haplotigs of one strain. Without aggressive
    //       cleaning, this information would be hard to use for the end user.
    //       (The other way around however, is easier, i.e. retrieving those
    //        variants from an aggressively cleaned long contig is relatively 
    //        trivial.)
    // HOW
    //     
    // TODO
    //     Linear buffer uses the simple stack right now and is ugly.
    // RETURN
    //     number of spots treated
    int ret = 0, nb_treated=0;
    int verbose = 0;

    int threshold_max_visit = 100;
    int threshold_max_utg_length = 100000;
    asg_t *auxsg = ug->g;
    uint32_t nv, wu;
    asg_arc_t *av;
    int flag;
    

    // init buffers
    stacku32_t ms;  // main stack
    stacku32_init(&ms);
    stacku32_t ms_block;  // main stack block offset and lengths packed as: offset<<16 | length
    stacku32_init(&ms_block);
    uint16_t block_offset, block_l, block_l_new;
    int idx;
    uint32_t handle=0, packed=0, sink=0;
    uint32_t *targets = (uint32_t*)malloc(sizeof(uint32_t) * threshold_max_visit);
    uint32_t *children = (uint32_t*)malloc(sizeof(uint32_t) * threshold_max_visit);
    int targets_n, children_n;

    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (ug->u.a[vu>>1].len<threshold_min_handle_length) continue;
        if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, -1)<=1) continue;
        
        // initial topo check: handle's immediate target shall have no backward branching
        if (verbose) {fprintf(stderr, "[debug::%s] handle %.6d, dir %d\n", __func__, (int)(vu>>1)+1, (int)(vu&1));}
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        flag = 1;
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) continue;
            if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, -1)>1) {
                flag = 0; 
                if (verbose) {fprintf(stderr, "[debug::%s]     violation: %.6d\n", __func__, (int)(av[i].v>>1)+1);}
                break;
            }
        }
        if (!flag) continue;

        // reset buffers
        stacku32_reset(&ms);
        stacku32_reset(&ms_block);
        stacku32_push(&ms, vu);
        stacku32_push(&ms_block, 1);

        flag = 1;
        // 1st pass: walk to see if we potentailly have a pop-able tangle 
        if (verbose) {fprintf(stderr, "[debug::%s]      1st pass\n", __func__);}
        while (ms.n<threshold_max_visit){
            stacku32_peek_last_item(&ms_block, &packed);
            block_offset = packed>>16;
            block_l = (uint16_t)packed;
            if (block_l==0){
                flag = 0; 
                if (verbose) {fprintf(stderr, "[debug::%s]      violation - empty block\n", __func__);}
                break;
            }
            block_l_new = 0;

            for (uint16_t i=block_offset; i<block_offset+block_l; i++){
                nv = asg_arc_n(auxsg, ms.a[i]);
                av = asg_arc_a(auxsg, ms.a[i]);
                for (int vui=0; vui<nv; vui++){
                    if (av[vui].del) continue;
                    if ( av[vui].v==vu ){  // shall never loop back to the handle
                        flag=0; 
                        if (verbose) {fprintf(stderr, "[debug::%s]      violation - loops to handle\n", __func__);}
                        break;
                    }  
                    if (!stack32_is_in_stack(&ms, av[vui].v)){
                        stacku32_push(&ms, av[vui].v);
                        block_l_new++;
                    }
                }
            }

            // check terminate: is there a unitig in the previous block that targets all unitigs added in the current batch
            if (ms_block.n>1){  // ==1 is the very start, don't need to check.
                if (verbose) {fprintf(stderr, "[debug::%s]      check termination (ms_block.n is %d)\n", __func__, ms_block.n);}
                // collect targets, it's a subset of the current batch 
                //  (if a unitig can be reached by another unitig in the current batch,
                //   we exclude the later.)
                int passed, yes_terminate;
                targets_n = 0;
                // fprintf(stderr, "sancheck: ms.n %d, ms_block.n %d; block_offset %d, block_l %d, block_l_new %d\n",
                //                 ms.n, ms_block.n, block_offset, block_l, block_l_new);
                for (uint16_t i=block_offset+block_l; i<block_offset+block_l+block_l_new; i++){
                    handle = ms.a[i];
                    passed = 1;
                    for (uint16_t j=block_offset+block_l; j<block_offset+block_l+block_l_new; j++){
                        if (j==i) continue;
                        if (hamt_check_if_is_immediate_decedent_of(auxsg, handle, ms.a[j])){
                            passed = 0;  // handle reaches another unitig in the current batch, ditch the handle.
                            break;
                        }
                    }
                    if (passed){
                        targets[targets_n] = handle;
                        targets_n++;
                    }
                }
                if (targets_n==0){
                    fprintf(stderr, "[W::%s] weird stuff\n", __func__);
                    flag = 0;
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s]        >targets: %d\n", __func__, targets_n);}
                    yes_terminate = 0;
                    for (uint16_t i=block_offset; i<block_offset+block_l; i++){
                        nv = asg_arc_n(auxsg, ms.a[i]);
                        av = asg_arc_a(auxsg, ms.a[i]);
                        children_n = 0;
                        for (uint32_t tmp=0; tmp<nv; tmp++){
                            if (av[tmp].del) continue;
                            children[children_n] = av[tmp].v;
                            children_n++;
                        }
                        if (uint32_buffer_unordered_equal(targets, targets_n, children, children_n)){
                            // ok terminate and go on to check more
                            sink = ms.a[i];
                            packed = ( ((uint32_t)(block_offset+block_l_new))<<16) | block_l_new;
                            stacku32_push(&ms_block, packed);
                            yes_terminate = 1;
                            if (verbose) {fprintf(stderr, "[debug::%s]      terminate, sink is %.6d\n", __func__, (int)(sink>>1)+1);}
                            break;
                        }   
                    }
                    if (yes_terminate) break;  // break out of 1st pass
                }
            }
            packed = ( ((uint32_t)(block_offset+block_l))<<16) | block_l_new;
            stacku32_push(&ms_block, packed);

            if (!flag) break;
        }
        if (verbose) {fprintf(stderr, "[debug::%s]      ms size %d ms_block size %d\n", __func__, ms.n, ms_block.n);}
        if (!flag) continue;

        // 2nd pass: check bifurcations
        if (verbose) {
            fprintf(stderr, "[debug::%s]      2nd pass, in the buffer: \n", __func__);
            for (int i=0; i<ms.n; i++){
                fprintf(stderr, "[debug::%s]        utg%.6d\n", __func__, (int)(ms.a[i]>>1)+1);
            }
        }
        flag = 1;
        for (int i=1; i<ms.n; i++){  // (skip the handle of the tangle)
            nv = asg_arc_n(auxsg, ms.a[i]^1);
            av = asg_arc_a(auxsg, ms.a[i]^1);
            for (int j=0; j<nv; j++){
                if (av[j].del) continue;
                if (av[j].v==(vu^1)) continue;
                if (!stack32_is_in_stack(&ms, av[j].v^1)){
                        if (verbose) {fprintf(stderr, "[debug::%s]      violation - utg %.6d backward bifurcation\n", __func__, (int)(av[j].v>>1)+1);}
                        flag = 0;
                        break;
                    }
            }
            // also, require that unitigs are short unless a) it's a tip, or b) it's the last batch.
            if (i<ms.n-1){
                if (ug->u.a[ms.a[i]>>1].len>threshold_max_utg_length && hamt_asgarc_util_countSuc(auxsg, ms.a[i], 0, 0, -1)!=0){
                    flag = 0;
                    if (verbose) {fprintf(stderr, "[debug::%s]      violation - touched a long non-tip unitig\n", __func__);}
                    break;
                }
            }
            if (!flag) break;
        }
        if (!flag) continue;


        // 3rd pass: get a path, drop everything else
        fprintf(stderr, "[debug::%s]      3rd pass i.e. going to pop tangle, sink is %.6d\n", __func__, (int)(sink>>1)+1);
        nb_treated += hamt_ug_pop_subgraph(sg, ug, sink^1, vu^1, 0, 1);
        ret+=1;

    }
    

    // free buffers
    stacku32_destroy(&ms);
    stacku32_destroy(&ms_block);
    free(targets);
    free(children);

    if (nb_treated){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    fprintf(stderr, "[M::%s] treated %d spots (%d unitigs)\n", __func__, ret, nb_treated);
    return ret;
}

// int hamt_ug_cutFalseLink_emulateForcDirectedLayout(){
//     // FUNC
//     //     The Raven assembler uses force-directed graph layout to find
//     //      enlarged arcs i.e. arcs that connect two otherwise topologically 
//     //      distant spots.
//     //     For sake of implementation simplicity, here we go after the topology 
//     //      without doing graph layout.

// }
int hamt_ug_popLooseTangles_v2(asg_t *sg, ma_ug_t *ug, int max_step){
    // FUNC
    //      Improved hamt_ug_popLooseTangles, no longer need the steps to be "synced".
    //      The flow:
    //        1) mark "strongly bonded components" for the whole graph; each component has at least 2 ends
    //           which are specially marked or stored,
    //        2) for each component try to pop through.
    // TODO
    //     Linear buffer uses the simple stack right now and is ugly.
    // RETURN
    //     number of spots treated
    int ret = 0, nb_treated=0;
    int verbose = 1;

    asg_t *auxsg = ug->g;
    uint32_t nv, vu, wu;
    asg_arc_t *av;
    
    int32_t *marks = (int32_t*)calloc(auxsg->n_seq, sizeof(int32_t));  // 0 is unset or doesn't matter; note that each unitig shall only have 1 label. This is enforce below, when we traverse the graph.
    uint8_t *specials = (uint8_t*)calloc(auxsg->n_seq, sizeof(uint8_t));  // 0 is unset, 1 is NOT special, 2 is special
    uint8_t is_special=2, is_not_special=1, is_unset=0;
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq*2, sizeof(uint8_t));  // for graph traversing
    uint8_t *node_is_weak = (uint8_t*)calloc(auxsg->n_seq, sizeof(uint8_t));  // if a node is not strongly bonded with all its neighbors (given no forbidden), then it's a weak node
    int32_t m = 1;
    int flag = 0;

    // step1a: mark weak nodes
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        // skip simple topo
        if ( (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0 && hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0) ||  // isolated unitig
             (hamt_asgarc_util_isIsolatedCircle(auxsg, vu, 0))  // isolated circle
             ){ 
            continue;
        }

        flag = 1;
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del) continue;
            if (!hamt_ug_check_if_arc_strongly_bonded(ug, vu, av[i].v, 0, max_step) || 
                !hamt_ug_check_if_arc_strongly_bonded(ug,av[i].v^1, vu^1, 0, max_step)){
                    node_is_weak[vu>>1] = 2;
                    if (verbose && flag) {fprintf(stderr, "[debug::%s] marked %.6d as weak node\n", __func__, (int)(vu>>1)+1);}
                    flag = 0;
                    node_is_weak[av[i].v>>1] = 2;
                    if (verbose) {fprintf(stderr, "[debug::%s] marked %.6d as weak node (via checking%.6d )\n", 
                                                    __func__, (int)(av[i].v>>1)+1,  (int)(vu>>1)+1);}   
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s] passed %.6d (checking against %.6d )\n", 
                                                    __func__, (int)(av[i].v>>1)+1,  (int)(vu>>1)+1);}   
            }
        }
    }
    // (if a node has weak nodes on both sides, mark it as weak too)
    while (1){
        flag = 0;
        int tmp = 0;
        for (vu=0; vu<auxsg->n_seq*2; vu++){
            if (node_is_weak[vu>>1]) continue;
            tmp = 0;
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if (node_is_weak[av[i].v>>1]>1){
                    tmp++;
                    break;
                }
            }
            if (!tmp) continue;

            flag++;
            node_is_weak[vu>>1] = 1;

            tmp = 0;
            nv = asg_arc_n(auxsg, vu^1);
            av = asg_arc_a(auxsg, vu^1);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if (node_is_weak[av[i].v>>1]>1){
                    tmp++;
                    break;
                }
            }
            if (!tmp) continue;

            node_is_weak[vu>>1] = 2;

            // flag++;
        }
        fprintf(stderr, "[debug::%s] flipped %d\n", __func__, flag);
        if (!flag) break;  // no more updates
    }
    for (vu=0; vu<auxsg->n_seq*2; vu++){  // debug, label by coloration in bandage and break out
        ug->utg_coverage[vu>>1] = (int)node_is_weak[vu>>1]*10 + 5;
    }
    return 1;



    ///////////////// step1b: mark components
    //         if any node of a strong arc is a weak node (via checking node_is_weak), the arc is considered to be weak .
    stacku32_t stack;  // DFS; BFS is also fine.
    stacku32_init(&stack);
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) continue;  // we will manually check both direction
        if (marks[vu>>1]!=0) continue;
        
        // early terminations
        if ( (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0 && hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0) ||  // isolated unitig
             (hamt_asgarc_util_isIsolatedCircle(auxsg, vu, 0))  // isolated circle
             ){ 
            marks[vu>>1] = m;
            m++;
            continue;
        }
        // skip, we will come to these from other places (or not, it doesn't matter)
        if (hamt_asgarc_util_isTip(auxsg, vu, 0, 0, 0)) continue;

        if (verbose){
            fprintf(stderr, "[debug::%s] stage1b, seed %.6d\n", __func__, (int)(vu>>1)+1);
        }
        // traverse
        memset(color, 0, auxsg->n_seq);
        stacku32_reset(&stack);
        stacku32_push(&stack, vu);
        stacku32_push(&stack, vu^1);
        while (stacku32_pop(&stack, &vu)){
            // if (marks[vu>>1]) {
                // assert(marks[vu>>1]==m);
                // fprintf(stderr, "[WARN] updated utg%.6d mark from %d to %d\n", (int)(vu>>1)+1, 
                //                             (int)marks[vu>>1], (int)m);
            // }
            marks[vu>>1] = m;
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (uint32_t i=0; i<nv; i++){
                if (av[i].del) continue;
                if (color[av[i].v]==0 && marks[av[i].v>>1]==0){
                    if (verbose) fprintf(stderr, "[debug::%s]        new target: %.6d\n", __func__, (int)(av[i].v>>1)+1);
                    if (hamt_ug_check_if_arc_strongly_bonded_bothdir(ug, vu, av[i].v, node_is_weak, max_step)){
                        // if (marks[av[i].v>>1]) {
                            // assert(marks[av[i].v>>1]==m);
                            // fprintf(stderr, "[WARN] updated utg%.6d mark from %d to %d\n", (int)(av[i].v>>1)+1, 
                            //                 (int)marks[av[i].v>>1], (int)m);
                        // }
                        marks[av[i].v>>1] = m;
                        stacku32_push(&stack, av[i].v);  // only push when the target is strongly bonded
                    }
                    color[av[i].v] = 1;  // mark a new discovery as seen, regardless of the bonding
                }
            }
            color[vu] = 2;
        }

        // update label
        m++;
    }

    // ~debug~
    // for (vu=0; vu<auxsg->n_seq*2; vu++){
    //     if (vu&1) continue;
    //     fprintf(stderr, "[strong bond mark] utg%.6d\t%d\n", (int)(vu>>1)+1, marks[vu>>1]);
    // }
    for (vu=0; vu<auxsg->n_seq*2; vu++){  // debug, label by coloration in bandage and break out
        ug->utg_coverage[vu>>1] = (int)marks[vu>>1];
    }
    return 1;

    // step2: mark special

    // step3: check each component and try popping

finish:
    // free buffers
    free(marks);
    free(specials);
    free(color);
    free(node_is_weak);
    stacku32_destroy(&stack);
    if (nb_treated){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    fprintf(stderr, "[M::%s] treated %d spots (%d unitigs)\n", __func__, ret, nb_treated);
    return ret;
}