#include <sys/stat.h>
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
#include "kvec.h"
#include "htab.h"
#include "kthread.h"
#include "t-sne.h"
#include "khashl.h"

#define HAMT_PI 3.141592653589793
#define HAMT_TWOPI 6.283185307179586
#define HAMT_MAX(x, y) ((x >= y)?(x):(y)) 
#define HAMT_MIN(x, y) ((x <= y)?(x):(y))
#define HAMT_DIFF(x, y) ((HAMT_MAX((x), (y))) - (HAMT_MIN((x), (y))))
#define UTG_LEN(ug, vu) ((ug)->u.a[(vu)>>1].len)  // vu has the direction bit

#define paf_ht_eq(a, b) ((a)==(b))
#define paf_ht_hash(a) ((a))
KHASHL_MAP_INIT(static klib_unused, paf_ht_t, paf_ht, uint64_t, uint16_t, paf_ht_hash, paf_ht_eq)
#define PAF_CT_PRE 12
#define PAF_CT_PRE_M ((1<<(PAF_CT_PRE))-1)
#define paf_ct_eq(a, b) ((a>>PAF_CT_PRE)==(b>>PAF_CT_PRE))
#define paf_ct_hash(a) ((a>>PAF_CT_PRE))
KHASHL_MAP_INIT(static klib_unused, paf_ct_t, paf_ct, uint64_t, uint32_t, paf_ct_hash, paf_ct_eq)

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)
KRADIX_SORT_INIT(ovhamt64, uint64_t, uint64_t, 8)
KRADIX_SORT_INIT(ovhamt32, uint32_t, uint32_t, 4)

typedef kvec_t(char*) kvec_strings_t;

#define HAMT_PRIMARY_LABEL 0
#define HAMT_ALTER_LABEL 1
void hamt_ug_util_BFS_markSubgraph_trailing(ma_ug_t *ug_old, ma_ug_t *ug_new, int base_label);

uint16_t index5NF[1024]={
0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 31, 
47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 15, 
62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 58, 
77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 43, 
92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 88, 103, 104, 105, 27, 
106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 73, 117, 118, 119, 11, 
120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 113, 131, 132, 133, 54, 
134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 99, 145, 146, 147, 39, 
148, 149, 150, 151, 152, 153, 154, 141, 155, 156, 157, 84, 158, 159, 160, 23, 
161, 162, 163, 164, 165, 166, 167, 127, 168, 169, 170, 69, 171, 172, 173, 7, 
174, 175, 176, 177, 178, 179, 180, 164, 181, 182, 183, 109, 184, 185, 186, 50, 
187, 188, 189, 190, 191, 192, 193, 151, 194, 195, 196, 95, 197, 198, 199, 35, 
200, 201, 202, 190, 203, 204, 205, 137, 206, 207, 208, 80, 209, 210, 211, 19, 
212, 213, 214, 177, 215, 216, 217, 123, 218, 219, 220, 65, 221, 222, 223, 3, 
224, 225, 226, 223, 227, 228, 229, 173, 230, 231, 232, 119, 233, 234, 235, 61, 
236, 237, 238, 211, 239, 240, 241, 160, 242, 243, 244, 105, 245, 246, 247, 46, 
248, 249, 250, 199, 251, 252, 253, 147, 254, 255, 256, 91, 257, 258, 247, 30, 
259, 260, 261, 186, 262, 263, 264, 133, 265, 266, 267, 76, 268, 269, 235, 14, 
270, 271, 272, 220, 273, 274, 275, 170, 276, 277, 278, 116, 279, 280, 267, 57, 
281, 282, 283, 208, 284, 285, 286, 157, 287, 288, 289, 102, 290, 291, 256, 42, 
292, 293, 294, 196, 295, 296, 297, 144, 298, 299, 289, 87, 300, 301, 244, 26, 
302, 303, 304, 183, 305, 306, 307, 130, 308, 309, 278, 72, 310, 311, 232, 10, 
312, 313, 314, 217, 315, 316, 317, 167, 318, 319, 307, 112, 320, 321, 264, 53, 
322, 323, 324, 205, 325, 326, 327, 154, 328, 329, 297, 98, 330, 331, 253, 38, 
332, 333, 334, 193, 335, 336, 327, 140, 337, 338, 286, 83, 339, 340, 241, 22, 
341, 342, 343, 180, 344, 345, 317, 126, 346, 347, 275, 68, 348, 349, 229, 6, 
350, 351, 352, 214, 353, 354, 343, 163, 355, 356, 304, 108, 357, 358, 261, 49, 
359, 360, 361, 202, 362, 363, 334, 150, 364, 365, 294, 94, 366, 367, 250, 34, 
368, 369, 361, 189, 370, 371, 324, 136, 372, 373, 283, 79, 374, 375, 238, 18, 
376, 377, 352, 176, 378, 379, 314, 122, 380, 381, 272, 64, 382, 383, 226, 2, 
384, 385, 383, 222, 386, 387, 349, 172, 388, 389, 311, 118, 390, 391, 269, 60, 
392, 393, 375, 210, 394, 395, 340, 159, 396, 397, 301, 104, 398, 399, 258, 45, 
400, 401, 367, 198, 402, 403, 331, 146, 404, 405, 291, 90, 406, 399, 246, 29, 
407, 408, 358, 185, 409, 410, 321, 132, 411, 412, 280, 75, 413, 391, 234, 13, 
414, 415, 381, 219, 416, 417, 347, 169, 418, 419, 309, 115, 420, 412, 266, 56, 
421, 422, 373, 207, 423, 424, 338, 156, 425, 426, 299, 101, 427, 405, 255, 41, 
428, 429, 365, 195, 430, 431, 329, 143, 432, 426, 288, 86, 433, 397, 243, 25, 
434, 435, 356, 182, 436, 437, 319, 129, 438, 419, 277, 71, 439, 389, 231, 9, 
440, 441, 379, 216, 442, 443, 345, 166, 444, 437, 306, 111, 445, 410, 263, 52, 
446, 447, 371, 204, 448, 449, 336, 153, 450, 431, 296, 97, 451, 403, 252, 37, 
452, 453, 363, 192, 454, 449, 326, 139, 455, 424, 285, 82, 456, 395, 240, 21, 
457, 458, 354, 179, 459, 443, 316, 125, 460, 417, 274, 67, 461, 387, 228, 5, 
462, 463, 377, 213, 464, 458, 342, 162, 465, 435, 303, 107, 466, 408, 260, 48, 
467, 468, 369, 201, 469, 453, 333, 149, 470, 429, 293, 93, 471, 401, 249, 33, 
472, 468, 360, 188, 473, 447, 323, 135, 474, 422, 282, 78, 475, 393, 237, 17, 
476, 463, 351, 175, 477, 441, 313, 121, 478, 415, 271, 63, 479, 385, 225, 1, 
480, 479, 382, 221, 481, 461, 348, 171, 482, 439, 310, 117, 483, 413, 268, 59, 
484, 475, 374, 209, 485, 456, 339, 158, 486, 433, 300, 103, 487, 406, 257, 44, 
488, 471, 366, 197, 489, 451, 330, 145, 490, 427, 290, 89, 487, 398, 245, 28, 
491, 466, 357, 184, 492, 445, 320, 131, 493, 420, 279, 74, 483, 390, 233, 12, 
494, 478, 380, 218, 495, 460, 346, 168, 496, 438, 308, 114, 493, 411, 265, 55, 
497, 474, 372, 206, 498, 455, 337, 155, 499, 432, 298, 100, 490, 404, 254, 40, 
500, 470, 364, 194, 501, 450, 328, 142, 499, 425, 287, 85, 486, 396, 242, 24, 
502, 465, 355, 181, 503, 444, 318, 128, 496, 418, 276, 70, 482, 388, 230, 8, 
504, 477, 378, 215, 505, 459, 344, 165, 503, 436, 305, 110, 492, 409, 262, 51, 
506, 473, 370, 203, 507, 454, 335, 152, 501, 430, 295, 96, 489, 402, 251, 36, 
508, 469, 362, 191, 507, 448, 325, 138, 498, 423, 284, 81, 485, 394, 239, 20, 
509, 464, 353, 178, 505, 442, 315, 124, 495, 416, 273, 66, 481, 386, 227, 4, 
510, 476, 376, 212, 509, 457, 341, 161, 502, 434, 302, 106, 491, 407, 259, 47, 
511, 472, 368, 200, 508, 452, 332, 148, 500, 428, 292, 92, 488, 400, 248, 32, 
511, 467, 359, 187, 506, 446, 322, 134, 497, 421, 281, 77, 484, 392, 236, 16, 
510, 462, 350, 174, 504, 440, 312, 120, 494, 414, 270, 62, 480, 384, 224, 0, 
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

char seqcmp(char n){
    if (n=='A') return 'T';
    if (n=='T') return 'A';
    if (n=='C') return 'G';
    if (n=='G') return 'C';
    return n;
}

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
int stacku32_peek_prev_item(stacku32_t *stack, int i, uint32_t *buf){
    // RETURN
    //     1 if success
    //     0 if invalid i (larger than stack size) or fail (stack is empty)
    if (stack->n>0 && i<stack->n){
        *buf = stack->a[i-1];
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



/**
 * @brief fixed length single bit array
*/
typedef struct {
    uint8_t *a;  // left to right
    uint32_t m, mm; // array size and expected bit-wise array size
}hamt_ba_t;

/**
 * @par bitsize the precise array length. Allocation is 8-bit blocks, 
 *      length sancheck is precise.
 * @func Allocate memory; will zero out.
*/
hamt_ba_t *hamt_ba_t_init(uint32_t bitsize){
    hamt_ba_t *d = (hamt_ba_t*)calloc(1, sizeof *d);
    uint32_t size = bitsize/8+1;
    d->a = (uint8_t*)calloc(size, 1);
    d->m = size;
    d->mm = bitsize;
    return d;
}
void hamt_ba_t_destroy(hamt_ba_t *d){
    free(d->a);
    free(d);
}

/**
 * @func Reallocate, will zero out new allocation. 
*/
void hamt_ba_t_resize(hamt_ba_t *d, uint32_t bitsize){
    uint32_t prev_size = d->m;
    uint8_t prev_last = d->a[d->m-1];
    d->m = bitsize/8 +1 ;
    d->mm = bitsize;
    d->a = (uint8_t*)realloc(d->a, d->m);
    if (!d->a){
        fprintf(stderr, "[E::%s] realloc failed\n", __func__);
        exit(1);
    }
    d->a[prev_size-1] = prev_last;
    memset(d->a+prev_size, 0, d->m - prev_size);
}

/**
 * @par i index in the bitarray
*/
uint8_t hamt_ba_t_get(hamt_ba_t *d, uint32_t i){
    if (i>d->mm){
        fprintf(stderr, "[E::%s] index exceeds array size\n", __func__);
        exit(1);
    }
    int idx = i/8;
    int shift = i%8;
    uint8_t ret = d->a[idx] & ((((uint8_t)1)<<7) >> shift);
    return ret>>(7-shift);
}

/**
 @par op 0 for unmark, 1 for mark, 2 for toggle
 @par i index in the bitarray
*/
void hamt_ba_t_write(hamt_ba_t *d, int op, uint32_t i){
    if (i>d->mm){
        fprintf(stderr, "[E::%s] index exceeds array size\n", __func__);
        exit(1);
    }
    int shift = i%8;
    int index = i/8;
    uint8_t mask = (((uint8_t)1)<<7) >> shift;
    if (op==0){
        d->a[index] &= (~mask);
    }else if (op==1){
        d->a[index] |= mask;
    }else if (op==2){
        d->a[index] ^= mask;
    }
}

/**
 * @brief print out bits
*/
void hamt_bat_t_debugprint(hamt_ba_t *d){
    for (int i=0; i<=d->mm; i++){
        fprintf(stderr, "[debug::%s] %d => %d\n", __func__, 
                                            i, (int)hamt_ba_t_get(d, i));
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

static inline void hamt_ug_remove_a_unitig(asg_t *sg, ma_ug_t *ug, uint32_t vu, int put_to_alt){
    // given vu, drop all of its arcs and move it to the alt graph (optional if put_to_alt)
    asg_t *auxsg = ug->g;
    uint32_t nv = asg_arc_n(auxsg, vu);
    asg_arc_t *av = asg_arc_a(auxsg, vu);
    for (int i=0; i<nv; i++){
        if (!av[i].del) hamt_ug_arc_del(sg, ug, vu, av[i].v, 1);
    }
    nv = asg_arc_n(auxsg, vu^1);
    av = asg_arc_a(auxsg, vu^1);
    for (int i=0; i<nv; i++){
        if (!av[i].del) hamt_ug_arc_del(sg, ug, vu^1, av[i].v, 1);
    }
    if (put_to_alt){
        hamt_ug_utg_softdel(sg, ug, vu, 1);
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
    stacku32_destroy(&stack);
    stacku64_destroy(&finishing_time);
}

/**
 * @brief Given the graph and two nodes (ordered), fetch the point to their edge.
 * @par resultp Output, is the address of point to be modified. 
 * @return 1 if found, 0 not found
*/
int hamt_get_arc(asg_t *g, uint32_t src, uint32_t dest, asg_arc_t** resultp){
    uint32_t i;
    if(g->seq[src>>1].del) return 0;

    uint32_t nv = asg_arc_n(g, src);
    asg_arc_t *av = asg_arc_a(g, src);
    for (i = 0; i < nv; i++)
    {
        if(av[i].del) continue;
        if(av[i].v == dest)
        {
            (*resultp) = &av[i];
            return 1;
        } 
    }

    return 0;
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
    //      Debug version of hamt_check_diploid.
    //      The denominator of ratio1 and ratio2 is min(len(vu1), len(vu2)).
    //       (see also: `hamt_check_diploid_report2`)
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
void hamt_check_diploid_report2(ma_ug_t *ug, uint32_t vu1, uint32_t vu2, ma_hit_t_alloc *reverse_sources,
                              float *ratio1, float *ratio2)
{
    // FUNC
    //      Debug version of hamt_check_diploid.
    //      The denominator of ratio1 and ratio2 is: len(vu1).
    //       (use `hamt_check_diploid_report` to use min(len(vu1), len(vu2))).
    int verbose = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu_short, vu_long;
    uint32_t qn, qn2, tn;
    int nb_targets;
    int found, cnt_hit=0;
    vecu32_t v;
    vecu32_init(&v);

    for (int i=0; i<ug->u.a[vu1>>1].n; i++){  // iterate over reads in the first query unitig
        qn = ug->u.a[vu1>>1].a[i]>>33;

        // check if there's any read in the other unitig that targets the read 
        found = 0;
        for (int j=0; j<ug->u.a[vu2>>1].n; j++){  // iterate over reads in the other unitig
            qn2 = ug->u.a[vu2>>1].a[j]>>33;
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
    
    *ratio1 = (float)cnt_hit / ug->u.a[vu1>>1].n;
    *ratio2 = (float)v.n / ug->u.a[vu1>>1].n;
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
    int verbose = 0;
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

/**
 * @brief Check if arc v->v exists, and v has no other arcs.
 * @return 1 if yes, 0 if no
*/
int hamt_asgarc_util_isSelfCirc(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
    if (hamt_asgarc_util_countSuc(g, v, include_del_seq, include_del_arc, base_label)==1){
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (int i=0; i<(int)nv; i++){
            if (!include_del_arc && av[i].del) {continue;}
            if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
            if (base_label>=0 && g->seq_vis[av[i].v>>1]!=base_label){continue;}
            if (av[i].v==v) return 1;
        }
    }
    return 0;
}


/**
 * @brief Check if v is a tip (either direction) with only 1 connection. 
 * @par v unitig ID, with the direction bit.
 * @return 1 if yes, 0 if no.
*/
int hamt_asgarc_util_isDanglingTip(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc, int base_label){
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

int hamt_asgarc_util_get_the_two_targets(asg_t *g, uint32_t v, uint32_t *w, uint32_t *u,
                                         int include_del_seq, int include_del_arc, int base_label){
    // NOTE
    //    Will not report error if there are more than two targets (will take the first two).
    // RET
    //    0 if ok
    //    1 if less than two targets
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
    if (idx<2){
        return 1;
    }else{
        *w = buf[0];
        *u = buf[1];
        return 0;
    }

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
    if (hamt_asgarc_util_get_the_two_targets(g, v2, &w2[0], &w2[1], 0, 0, base_label)!=0)
        return 0;
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

/**
 * @return 1 if yes, 0 if no
*/
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
    if (verbose){
        fprintf(stderr, "[debug::%s] utg%.6d (%d) vs utg%.6d (%d)\n", __func__, 
            (int)(v0>>1)+1, (int)(v0&1), (int)(w0>>1)+1, (int)(w0&1));
    }

    // if ((v0>>1)+1==993 || (v0>>1)+1==95) verbose = 2;  // 990

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
        if (ug->u.a[source_nodir].len<100000){
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

        if ((max_visited>0) && (nb_visited>max_visited) ){  // search span too big, assume failure
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
        fprintf(stderr, "[E::%s]ERROR u1\n", __func__);
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
            fprintf(stderr, "[E::%s]error1 \n", __func__);
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
                    fprintf(stderr, "[E::%s]ERROR, tmp_v1\n", __func__);
                    exit(1);
                }
                if (hamt_asgarc_util_get_the_one_target(auxsg, w2[1], &tmp_v2, 0, 0, base_label)<0){
                    fprintf(stderr, "[E::%s]ERROR tmp_v2\n", __func__);
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
        if (hamt_asgarc_util_get_the_two_targets(auxsg, v0, &v1, &v2, 0, 0, base_label)!=0)
            continue;
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
        if (hamt_asgarc_util_get_the_two_targets(auxsg, v3, &v_tmp1, &v_tmp2, 0, 0, base_label)!=0)
            continue;
        if ((v_tmp1>>1)==(v1>>1) || (v_tmp1>>1)==(v2>>1)){
            v4 = v_tmp2;
        }else{
            v4 = v_tmp1;
        }
        if (hamt_asgarc_util_get_the_two_targets(auxsg, v3^1, &v_tmp1, &v_tmp2, 0, 0, base_label)!=0)
            continue;
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
                fprintf(stderr, "[E::%s] ERROR\n", __func__);
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


/**
 * @brief Helper function. Treats the bidirectional graph as undirected, find
 *        disconnected subgraphs, and label their nodes with unique non-negative 
 *        integer IDs.
 * @par sg String graph or aux grpah of a unitig/contig graph. Note that for string graph,
 *        set base_label to -1.
 * @par max_label Output, stores the number of disconnected subgraphs.
 * @par base_label 0 for primary , 1 for alt, -1 ignores the check.
 * @return An integer array of length n_seq (no direction bit).
*/
int *hamt_asgarc_util_BFS_markSubgraph_core(asg_t *sg, int *n_subgraphs, int base_label){
    int verbose = 0;

    uint32_t v0, v, nv, w;
    asg_arc_t *av;
    int cnt_nread, cnt_bp;
    int tmp;
    int label = 0;
    int n_nodes;

    queue32_t q;
    queue32_init(&q);
    uint8_t *color = (uint8_t*)calloc(sg->n_seq, 1);
    int *subg_labels = (int*)calloc(sg->n_seq, sizeof(int));

    for (v0=0; v0<sg->n_seq*2; v0+=2){
        if (color[v0>>1]!=0) continue;
        if (base_label>=0 && sg->seq_vis[v0>>1]!=base_label){continue;}
        if (verbose){
            fprintf(stderr, "[debug::%s] seed utg%.6d\n", __func__, (int)(v0>>1)+1);
        }
        
        // BFS, treat graph as undirected
        n_nodes = 0;
        v = v0>>1;
        queue32_reset(&q);
        queue32_enqueue(&q, v);
        color[v] = 1;
        n_nodes++;
        subg_labels[v] = label;
        while(queue32_dequeue(&q, &v)){
            if (color[v]==2){
                continue;
            }
            
            for (uint32_t dir=0; dir<2; dir++){
                nv = asg_arc_n(sg, (v<<1)|dir);
                av = asg_arc_a(sg, (v<<1)|dir);
                for (uint32_t i=0; i<nv; i++){
                    if (av[i].del){continue;}
                    if (base_label>=0 && sg->seq_vis[av[i].v>>1]!=base_label) {continue;}
                    w = av[i].v>>1;
                    if (color[w]==0){
                        queue32_enqueue(&q, w);
                        color[w] = 1;
                        n_nodes++;
                        subg_labels[w] = label;
                    }
                }
                color[v] = 2;
            }
        }
        label++;
        if (verbose){
            fprintf(stderr, "[debug::%s]   got subgraph #%d, has %d nodes.\n", __func__, label, n_nodes);
        }
    }
    queue32_destroy(&q);
    free(color);
    *n_subgraphs = label;
    return subg_labels;
}

int hamt_ug_util_BFS_markSubgraph(ma_ug_t *ug, int base_label){
    // FUNC
    //     Traverse the ug and give each unitig a subgraph ID (arbitrary; 
    //       stable in terms of assembly runs, NOT stable between ug generations).
    // RET
    //     The largest subgraph ID.
    asg_t *auxsg = ug->g;
    int max_label;
    int *labels = hamt_asgarc_util_BFS_markSubgraph_core(auxsg, &max_label, -1);
    for (uint32_t i=0; i<auxsg->n_seq*2; i+=2){
        ug->u.a[i>>1].subg_label = labels[i>>1];
    }
    free(labels);
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

    g->arc = (asg_arc_t*)calloc(g0->n_arc, sizeof(asg_arc_t));
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
    uint32_t n_del=0, n_del_tot=0;
    int pre, suc, l, ret;

    asg64_v a = {0,0,0};
    int i_round;
    
    for (i_round=0; i_round<3; i_round++)
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
        }
        
        n_del += hamt_sg_pop_simpleInvertBubble(g);
        n_del_tot += n_del;
    }
    free(a.a);

    fprintf(stderr, "[M::%s] did %d rounds, dropped %d spots, used %.1f s\n\n", 
                        __func__, i_round, n_del_tot, Get_T()-startTime);

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

    int verbose = 0;
    // int minl = 50000;
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
        // if (ug->u.a[v>>1].len<minl || ug->u.a[w>>1].len<minl) continue;

        cov[2] = 0;
        av = asg_arc_a(auxsg, w^1);
        for (uint32_t i=0; i<asg_arc_n(auxsg, w^1); i++){  // collect the largest coverage of w's predecessor (except v which will be accounted by cov[0])
            if ((av[i].v>>1)==(v>>1)){  // ignore self
                continue;
            }
            if (av[i].del){continue;}
            if (base_label>=0 && (auxsg->seq_vis[av[i].v>>1]!=base_label)) {continue;}
            // if (ug->u.a[av[i].v>>1].len<minl) continue;

            covtmp = ug->utg_coverage[av[i].v>>1];
            if (covtmp>cov[2]){
                cov[2] = covtmp;
            }
        }
        if (cov[2]==0) continue;  // no child node were checked
                                  //
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
                if (verbose) fprintf(stderr, "[debug::%s] ignore\n", __func__);
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
    int ignore_by_length = 0, ignore_by_target = 0, spared = 0;

    asg_t *auxsg = ug->g;
    uint32_t vu, wu, nv;
    asg_arc_t *av;
    int nb_cut = 0;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (base_label>=0 && auxsg->seq_vis[vu>>1]!=base_label) {continue;}
        if (ug->u.a[vu>>1].len>100000){  // tip not short
            ignore_by_length++;
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
        for (uint32_t i=0; i<nv; i++){
            if (av[i].del){continue;}
            wu = av[i].v;
            // if (hamt_asgarc_util_countNoneTipSuc(auxsg, wu^1, base_label)==0){
            if (hamt_asgarc_util_countSuc(auxsg, wu^1, 0, 0, base_label)==1){
                if (verbose){
                    fprintf(stderr, "[debug::%s] spared an arc for tip utg%.6d (target utg%.6d)\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
                }
                spared++;
            }else{  // cut
                hamt_ug_arc_del(sg, ug, vu, wu, 1);
                nb_cut++;
                if (verbose){
                    fprintf(stderr, "[debug::%s] cut tip: utg%.6d \n", __func__, (int)(vu>>1)+1);
                }
            }
        }
        #endif


    }

    if (nb_cut){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    if (VERBOSE){
        fprintf(stderr, "[M::%s] cut %d unitig tigs; ignored not short %d , bifur %d, spared %d.\n", 
                            __func__, nb_cut, ignore_by_length, ignore_by_target, spared);
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
    if (hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label)!=0){
        if (verbose){
            fprintf(stderr, "[debug::%s] utg %.6d , failed to get wu0 and wu1\n", __func__, (int)(vu>>1)+1);
        }
        return 0;
    }
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



void hamt_asgarc_markBridges_DFS_v0(asg_t *sg, uint32_t v, uint32_t v_parent, int *DFStime,
                                 uint8_t *visited, int *tin, int *low, int base_label){
    // recursive.
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
            hamt_asgarc_markBridges_DFS_v0(sg, w, v, DFStime, visited, tin, low, base_label);
            low[v>>1] = low[v>>1]<low[w>>1]? low[v>>1] : low[w>>1];
            if (low[w>>1]>tin[v>>1]){  // arc is bridge
                av[i].is_bridge = 1;
                // mark the complementary arc as bridge too
                sg->arc[asg_arc_get_complementaryID(sg, &av[i])].is_bridge = 1;
            }
        }
    }
}

/**
 * @brief Treat the assembly graph as an undirected graph and mark bridge edges.
 *        Non-recursive version.
 * @par sg Aux graph of a unitig or contig graph. Will modify sg->arc[i].is_bridge.
 * @par root A node as DFS seed, has direction bit.
 * @par visited An array of length n_seq provided by caller.
 * @par tin An array of length n_seq initialized all to -1 provided by caller  
 * @par low (Similar to tin.)
 * @par base_label -1 to ignore, 0 to use primary graph, 1 to use alt graph.
 * @return Number of bridges.
*/
int hamt_asgarc_markBridges_DFS(asg_t *sg, uint32_t root,
                                uint8_t *visited, int *tin, int *low, int base_label) {
    int verbose = 0;
    int nb_bridge = 0;
    
    // for traversing graph
    int DFStime = 0;
    uint32_t vu, wu, parent, vu_;
    asg_arc_t *av_fwd, *av_bwd, *av;  
    uint32_t nv_fwd, nv_bwd;
    stacku32_t DFSstack;
    stacku32_init(&DFSstack);
    // for fetching arcs
    asg_arc_t *av0, *avr, *avtmp;  

    // init
    stacku32_push(&DFSstack, root);
    visited[root >> 1] = 1;
    tin[root >> 1] = DFStime;
    low[root >> 1] = tin[root >> 1];
    DFStime += 1;

    int exhausted;
    uint8_t parent_sancheck;
    int hasDFSparent;
    while (stacku32_pop(&DFSstack, &vu)) {
        if (verbose>1) {fprintf(stderr, "[debug::%s] pop vu tig %.6d\n", __func__, (int)(vu>>1)+1);}
        av_fwd = asg_arc_a(sg, vu);
        av_bwd = asg_arc_a(sg, vu ^ 1);
        nv_fwd = asg_arc_n(sg, vu);
        nv_bwd = asg_arc_n(sg, vu ^ 1);
        exhausted = 1;
        hasDFSparent = stacku32_peek_last_item(&DFSstack, &parent);
        for (uint32_t i_ = 0, i = 0; i_ < (nv_fwd + nv_bwd); i_++) {  // check all neighbors
            if (i_ < nv_fwd) {
                i = i_;
                av = av_fwd;
                vu_ = vu;
            } else {
                i = i_ - nv_fwd;
                av = av_bwd;
                vu_ = vu^1;
            }
            if (av[i].del) { continue; }
            if (base_label >= 0 && sg->seq_vis[av[i].v >> 1] != base_label) { continue; }

            wu = av[i].v;
            if (visited[wu >> 1] == 0) {  // found white
                if (verbose>1) {fprintf(stderr, "[debug::%s]   push wu tig %.6d\n", __func__, (int)(wu>>1)+1);}
                exhausted = 0;
                stacku32_push(&DFSstack, vu);  // put current node back
                stacku32_push(&DFSstack, wu);
                visited[wu>>1] = 1;
                DFStime++;
                low[wu >> 1] = DFStime;
                tin[wu >> 1] = DFStime;  // discovery time
                break;  // step deeper
            } else{
                if (hasDFSparent && (wu>>1)!=(parent>>1))
                    low[vu>>1] = low[vu>>1]<low[wu>>1]? low[vu>>1] : low[wu>>1];
            }
        }
        if (hasDFSparent==0) continue;
        if (exhausted && DFSstack.n!=0){  // update times
            if (verbose>1){
                fprintf(stderr, "[debug::%s] finish %.6d , DFSparent=%.6d, tin=%d low=%d\n", __func__,
                                (int)(vu>>1)+1, (int)(parent>>1)+1, 
                                tin[vu>>1], low[vu>>1]);
            }

            if (tin[vu >> 1] == low[vu >> 1]) {  // the edge between this node and its DFS parent is a bridge
                if (hamt_asgarc_util_isTip(sg, vu, 0, 0, base_label) ||  
                    hamt_asgarc_util_isTip(sg, parent, 0, 0, base_label)){
                    if (verbose>1) {fprintf(stderr, "[debug::%s] IGNORE because tip: edge %.6d to %.6d\n",
                                    __func__, (int)(vu>>1)+1, (int)(parent>>1)+1);}
                }else{
                    int found1 = hamt_get_arc(sg, parent,vu,  &av0);
                    if (!found1){
                        if (verbose>1) fprintf(stderr, "flip parent direction because found1 is false\n");
                        parent ^=1;
                        found1 = hamt_get_arc(sg, parent,vu,  &av0);
                    }
                    int found2 = hamt_get_arc(sg,vu, parent, &avtmp);  
                    int found3 = hamt_get_arc(sg,vu^1, parent^1, &avr);  
                    
                    if (verbose>1) {
                        fprintf(stderr, "vu %.6d parent %.6d stat1 %d stat2 %d stat3 %d\n", 
                                (int)(vu>>1)+1, (int)(parent>>1)+1, 
                                found1, found2, found3);
                        if (found1) fprintf(stderr, "av0: %.6d %.6d\n", (int)(av0->v>>1)+1, (int)(av0->ul>>33)+1);
                        if (found2) fprintf(stderr, "avtmp: %.6d %.6d\n", (int)(avtmp->v>>1)+1, (int)(avtmp->ul>>33)+1);
                        if (found3) fprintf(stderr, "avr: %.6d %.6d\n", (int)(avr->v>>1)+1, (int)(avr->ul>>33)+1);
                    }

                    if (found2) continue;

                    if (!found3){
                        fprintf(stderr, "[E::%s] failed to get the reverse arc: parent=%.6d (%d), vu=%.6d (%d)\n", 
                                __func__, (int)(parent>>1)+1, (int)(parent&1),
                                (int)(vu>>1)+1, (int)(vu&1));
                    }else{
                        av0->is_bridge = 1;
                        avr->is_bridge = 1;
                        nb_bridge++;
                    }
                    if (verbose) {
                        fprintf(stderr, "[debug::%s] MARK: edge %.6d (%d) to %.6d (%d)\n",
                                    __func__, (int)(vu>>1)+1, (int)(vu&1),
                                    (int)(parent>>1)+1, (int)(parent&1));
                    }
                }
            }
        }
    }
    stacku32_destroy(&DFSstack);
    return nb_bridge;
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
    int *tin = (int*)calloc(auxsg->n_seq, sizeof(int));
    int *low = (int*)calloc(auxsg->n_seq, sizeof(int));
    memset(tin, -1, auxsg->n_seq*sizeof(int));
    memset(low, -1, auxsg->n_seq*sizeof(int));
    int DFStime = 0;
    for (v=0; v<auxsg->n_seq*2; v++){
        if (v&1) continue;
        if (visited[v>>1]==0){
            // hamt_asgarc_markBridges_DFS(auxsg, v, v, &DFStime, visited, tin, low, base_label);
            hamt_asgarc_markBridges_DFS(auxsg, v, visited, tin, low, base_label);
            memset(tin, -1, auxsg->n_seq*sizeof(int));
            memset(low, -1, auxsg->n_seq*sizeof(int));
        }
    }

    // check arcs
    // (note: since the asg is treated as if undirected, edges in bubbles are not bridges)
    int nb_bridge = 0;
    for (uint32_t i=0; i<auxsg->n_arc; i++){
        if (auxsg->arc[i].del){continue;}
        if (!auxsg->arc[i].is_bridge){
            continue;
        }
        nb_bridge++;

        v = auxsg->arc[i].ul>>32;
        w = auxsg->arc[i].v;
        if (base_label>=0 && ( auxsg->seq_vis[v>>1]!=base_label || auxsg->seq_vis[w>>1]!=base_label) ) {continue;}
        if (verbose>1) {fprintf(stderr, "[debug::%s] checking at utg%.6d vs utg%.6d\n", __func__, (int)(v>>1)+1, (int)(w>>1)+1);}

        // check node length, no trust for coverage estimation in short nodes
        // if (ug->u.a[v>>1].len<50000 && ug->u.a[w>>1].len<50000) continue;

        // check node coverage
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
            if (hamt_asgarc_util_isTip(auxsg, v, 0, 0, base_label) &&
                (hamt_asgarc_util_countSuc(auxsg, v, 0, 0, base_label)==1)){
                // don't remove tips from circles
                // (however, allow cutting if the tip is still connected to somewhere else)
                if (verbose>1) {fprintf(stderr, "[debug::%s]    spared bcs handle is a tip\n", __func__);}
                continue;
            }

            // cut
            hamt_ug_arc_del(sg, ug, v, w, 1);  // deletes both direction
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

    fprintf(stderr, "[M::%s] cut %d bridges outof %d, used %.1fs.\n", 
            __func__, nb_cut, (nb_bridge-nb_cut)/2, Get_T()-startTime);
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
    //     Core function expects the path to be: start->...-> end (not the end^1 closing)
    //     Topo context of vu will NOT be checked; caller is responsible for sanchecks.
    //     Will search more than start/end vertex - since we can't add arc if the overlap was containment etc.
    //       if all recoverable overlaps were not ideal, do nothing.
    //     After recovering the arc, also trim off reads involved in containment. (reads will be DELETED)
    // RETURN
    //     1 upon recovering an arc (+cleanup)
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

    if (verbose) {fprintf(stderr, "[debug::%s] check vu utg%.6d dir %d\n", __func__, (int)(v0>>1)+1, (int)(v0&1));}

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
            break;
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
        stacku32_destroy(&DFSstack);
        stacku32_destroy(&alt_buf);
        free(packed);
        free(color);

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
    stacku32_destroy(&DFSstack);
    stacku32_destroy(&alt_buf);
    free(packed);
    free(scores); free(pis);
    free(color);

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
    fprintf(stderr, "[M::%s] done\n", __func__);
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


int hamt_ug_drop_redundant_nodes(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label){
    // FUNC
    //     We remove a node `v` if:
    //      1) there's at least one node `w` that has all predecessors of v
    //      2) also all sucessors of v
    //      3) and that all sucessors of v shall not have backward branchings that's not covered by remaining...(TODO)
    //     It is hard to correctly check 3) if there's a tangle. We just examine the simplest case.
    //     Implementation is handle-based, this function might need to be called multiple times.
    // RETURN
    //     # of nodes dropped.
    int total = 0;
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t nv, nw; 
    asg_arc_t *av, *aw;
    stacku32_t targets, will_be_covered;
    stacku32_init(&targets);
    stacku32_init(&will_be_covered);
    // pack target counts and unitig ID for sorting
    // note: maybe we can use uint32_t (since # unitigs is much smaller than # reads in regular datasets, 
    //       and # reads is less than 1<<28 instead of 32 in the current implementation).
    //       using u64 because it's simpler and either way it doesn't matter.
    stacku64_t buf_for_sort;  
    stacku64_init(&buf_for_sort); 

    int ret = 1, iter=0;
    while (ret && iter<50){
        ret = 0;
        for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
            if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)<=1) continue;
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            int i=0; 
            stacku32_reset(&targets);
            stacku32_reset(&will_be_covered);
            stacku64_reset(&buf_for_sort);
            for (i=0; i<nv; i++){
                if (av[i].del) continue;
                nw = asg_arc_n(auxsg, av[i].v);
                aw = asg_arc_a(auxsg, av[i].v);
                if (hamt_asgarc_util_countPre(auxsg, av[i].v, 0, 0, base_label)>1) {
                    // backward bifur, we won't drop this unitig, so all its targets will be reacheable.
                    for (int j=0; j<nw; j++){
                        if (aw[j].del) continue;
                        stacku32_push(&will_be_covered, aw[j].v);
                    }
                }else{
                    // collect targets. Don't treat tips differently.
                    for (int j=0; j<nw; j++){
                        if (aw[j].del) continue;
                        if (!stack32_is_in_stack(&targets, aw[j].v)) stacku32_push(&targets, aw[j].v);
                    }
                    // put the current unitig into popping queue 
                    stacku64_push(&buf_for_sort, (((uint64_t)(hamt_asgarc_util_countSuc(auxsg, av[i].v, 0, 0, base_label)))<<32) | av[i].v);
                }
            }
            // attach a del bit: some targets do not need to be covered, as their origins won't be popped.
            int left = 0;
            for (i=0; i<targets.n; i++){
                if (stack32_is_in_stack(&will_be_covered, targets.a[i])){
                    targets.a[i] = (targets.a[i]<<1) | 1;  // ok to no cover this entry
                }else{
                    targets.a[i] = (targets.a[i]<<1) | 0;
                    left++;
                }
            }
            if (buf_for_sort.n==0) continue;  // nothing to consider to be dropped
            
            // Sort. 
            radix_sort_ovhamt64(buf_for_sort.a, buf_for_sort.a+buf_for_sort.n);

            // Instead of solving a set covering problem, we simply check from the origin unitigs with (roughly) the most
            //  target counts, retain until all targets are covered.
            // "Roughly" because we do not care about some targets, but counting did not consider this. Doesn't matter too much.
            for (i=buf_for_sort.n-1; i>=0; i--){
                uint32_t wu = (uint32_t) buf_for_sort.a[i];
                if (left<=0){  // all covered, drop encountered origin unitig
                    if (left<0) {
                        fprintf(stderr, "[E::%s] minus leftover (continue anyway)\n", __func__);
                    }
                    if (ug->u.a[vu>>1].len>size_limit_bp) continue;
                    asg_arc_t *tmp = asg_arc_a(auxsg, wu);
                    for (int ii=0; ii<asg_arc_n(auxsg, wu); ii++){
                        hamt_ug_arc_del(sg, ug, wu, tmp[ii].v, 1);
                    }
                    tmp = asg_arc_a(auxsg, wu^1);
                    for (int ii=0; ii<asg_arc_n(auxsg, wu^1); ii++){
                        hamt_ug_arc_del(sg, ug, wu^1, tmp[ii].v, 1);
                    }
                    hamt_ug_utg_softdel(sg, ug, wu, 1);
                    ret++;
                    if (verbose) {
                        fprintf(stderr, "[debug::%s] handle utg%.6d drop %.6d\n", __func__, 
                                            (int)(vu>>1)+1, (int)(wu>>1)+1);
                    }
                }else{
                    nw = asg_arc_n(auxsg, wu);
                    aw = asg_arc_a(auxsg, wu);
                    for (int j=0; j<nw; j++){
                        if (aw[j].del) continue;
                        int idx = stacku32_index_value(&targets, aw[j].v<<1);
                        if (idx!=-1){
                            left--;
                            targets.a[idx] |=1;
                        }
                    }
                }
            }
        }
        total+=ret;
        iter++;
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        if (VERBOSE){fprintf(stderr, "[M::%s]    iter %d dropped %d\n", __func__, iter, ret);}
    }

    stacku32_destroy(&targets);
    stacku32_destroy(&will_be_covered);
    stacku64_destroy(&buf_for_sort);

    if (VERBOSE){fprintf(stderr, "[M::%s] dropped %d\n", __func__, total);}
    return ret;
}
int hamt_ug_drop_redundant_nodes_bruteforce(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label, int verbose){
    // FUNC
    //    Instead of checking locally, this function finds all "equivalent" or "sub-squivalent" unitig pairs.
    //    "Equivalent" means having identical set of predecessors and identical targets (directional).
    //    "Sub-equivalent" means unitigA has subsets of unitigB's predecessors and targets (directional),
    //      where we can drop A for almost no loss.
    // RET
    //    # of treatments
    int ret = 1, total = 0, iter=0;
    // int verbose = 0;
    float max_diff_ratio = 5;  // do no treat a pair if one of the unitig is significantly longer than the other one.

    asg_t *auxsg = ug->g;
    uint32_t nv;
    asg_arc_t *av;
    stacku32_t *h;
    stacku32_t candidates;
    stacku32_init(&candidates);

    stacku32_t *buf = (stacku32_t*)malloc(sizeof(stacku32_t) * auxsg->n_seq*2);
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){
        stacku32_init(&buf[i]);
    }


    while (ret && iter<50){
        ret = 0;
        for (uint32_t i=0; i<auxsg->n_seq*2; i++){
            stacku32_reset(&buf[i]);
        }   
        // Collect predecessors and targets. Do not treat tips differently when they are predecessors/targets.
        // "Packing" in the stack is:
        //       {0,1}      |  predecessors....  | targets... | nb_predecessors | nb_targets
        //  ^^1 if deleted                      
        //  (directions of the predecessors are: pre -> unitig -> targets.)
        for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
            h = &buf[vu];
            int cnt1=0, cnt2=0;

            // ignore tips or unconnected unitigs
            if (hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)==0 || 
                hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)==0){
                    stacku32_push(h, (uint32_t)1);
                    continue;  
            }

            // all green, collect info
            stacku32_push(h, (uint32_t)0);

            nv = asg_arc_n(auxsg, vu^1);
            av = asg_arc_a(auxsg, vu^1);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                stacku32_push(h, av[i].v^1);
                cnt1++;
            }

            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                stacku32_push(h, av[i].v);
                cnt2++;
            }
            
            stacku32_push(h, (uint32_t)cnt1);
            stacku32_push(h, (uint32_t)cnt2);
        }


        // Compare. When there's a hit that fufills all requirements, 
        //   we keep the one with higher coverage.
        // The requirements: (sub)equivalence, max length, max length diff ratio.
        for (uint32_t vu=0; vu<auxsg->n_seq*2-1; vu++){
            if (buf[vu].a[0]) continue;  // vu is ignored or has been removed
            // collect candidates            
            stacku32_reset(&candidates);
            //   (dir 1)
            nv = asg_arc_n(auxsg, vu^1);
            av = asg_arc_a(auxsg, vu^1);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if ((av[i].v>>1)==(vu>>1)) continue;
                int nw = asg_arc_n(auxsg, av[i].v^1);
                asg_arc_t *aw = asg_arc_a(auxsg, av[i].v^1);
                for (int j=0; j<nw; j++){
                    if (aw[j].del) continue;
                    if ((aw[j].v>>1)==(vu>>1)) continue;
                    if (!stack32_is_in_stack(&candidates, aw[j].v)) stacku32_push(&candidates, aw[j].v);
                }
            }
            //   (dir 2)
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (int i=0; i<nv; i++){
                if (av[i].del) continue;
                if ((av[i].v>>1)==(vu>>1)) continue;
                int nw = asg_arc_n(auxsg, av[i].v^1);
                asg_arc_t *aw = asg_arc_a(auxsg, av[i].v^1);
                for (int j=0; j<nw; j++){
                    if (aw[j].del) continue;
                    if ((aw[j].v>>1)==(vu>>1)) continue;
                    if (!stack32_is_in_stack(&candidates, aw[j].v^1)) stacku32_push(&candidates, aw[j].v^1);
                }
            }

            for (int i_wu=0; i_wu<candidates.n; i_wu++){
                uint32_t wu = candidates.a[i_wu];
                if (buf[wu].a[0]) continue;  // wu is not considered (e.g. it's a tig) or has been removed
                if ((vu>>1)==(wu>>1)) continue;   // skip, self comparison
                if (verbose){
                    fprintf(stderr, "[debug::%s] check utg%.6d (dir %d) and utg%.6d (dir %d)\n", 
                                __func__, (int)(vu>>1)+1, (int)(vu&1),
                                        (int)(wu>>1)+1, (int)(wu&1));
                }

                
                float len_diff_ratio = ((float)ug->u.a[vu>>1].len)/((float)ug->u.a[wu>>1].len);
                len_diff_ratio = len_diff_ratio<1? 1/len_diff_ratio : len_diff_ratio;
                if (len_diff_ratio>max_diff_ratio) continue;  // skip, length difference is too large
                if (verbose){fprintf(stderr, "[debug::%s]     length ratio ok\n", __func__);}

                if (ug->u.a[vu>>1].len>size_limit_bp || ug->u.a[wu>>1].len>size_limit_bp) continue;  // skip, one of the unitigs is too large
                if (verbose){fprintf(stderr, "[debug::%s]     max length ok\n", __func__);}

                int n_pre_vu = buf[vu].a[buf[vu].n-2];
                int n_suc_vu = buf[vu].a[buf[vu].n-1];
                int n_pre_wu = buf[wu].a[buf[wu].n-2];
                int n_suc_wu = buf[wu].a[buf[wu].n-1];
                int n_pre_min, n_suc_min;  // lengths of the smaller buffers
                int n_pre_max, n_suc_max;  // lengths of the larger buffers
                uint32_t the_larger, the_smaller;

                if (n_pre_vu>=n_pre_wu && n_suc_vu>=n_suc_wu){
                    n_pre_max = n_pre_vu;
                    n_suc_max = n_suc_vu;
                    n_pre_min = n_pre_wu;
                    n_suc_min = n_suc_wu;
                    the_larger = vu;
                    the_smaller = wu;
                }else if (n_pre_vu<=n_pre_wu && n_suc_vu<=n_suc_wu){
                    n_pre_max = n_pre_wu;
                    n_suc_max = n_suc_wu;
                    n_pre_min = n_pre_vu;
                    n_suc_min = n_suc_vu;
                    the_larger = wu;
                    the_smaller = vu;
                }else{  // can't be a hit
                    continue;
                }
                if (verbose){fprintf(stderr, "[debug::%s]     suc/tar sizes ok\n", __func__);}
                
                
                // expect all entries in the smaller buffer to be found in the other larger or equal sized buffer
                int passed = 1;
                for (int n=0; n<n_pre_min; n++){
                    if (!stack32_is_in_stack_givenrange(&buf[the_larger], buf[the_smaller].a[1+n], 1, 1+n_pre_max)){
                        passed = 0;
                        break;
                    }
                }
                if (!passed) continue;
                if (verbose){fprintf(stderr, "[debug::%s]     suc ok\n", __func__);}
                for (int n=0; n<n_suc_min; n++){
                    if (!stack32_is_in_stack_givenrange(&buf[the_larger], buf[the_smaller].a[1+n_pre_min+n], 1+n_pre_max, 1+n_pre_max+n_suc_max)){
                        passed = 0;
                        break;
                    }
                }
                if (!passed) continue;
                if (verbose){fprintf(stderr, "[debug::%s]     tar ok\n", __func__);}

                if (passed){  // try to drop one of the unitig 
                    uint32_t drop = ug->utg_coverage[vu>>1]<ug->utg_coverage[wu>>1]? vu : wu;
                    hamt_ug_remove_a_unitig(sg, ug, drop, 1);
                    buf[drop].a[0] = 1;
                    ret++;
                    if (verbose) {fprintf(stderr, "[debug::%s]     all green, drop utg%.6d\n", __func__, (int)(drop>>1)+1);}

                    // leave the inner loop if we decided to drop vu (which is the handle given by the outer loop)
                    // note/TODO: this is arbitrarily dependent on the order of checking. 
                    //            To do better, maybe we can check all pairs and log the decisions,
                    //             then find a set of cuts that maximizes the per base coverage after cutting.
                    if (drop==vu) break;  
                }
            }
        }

        iter++;
        total+=ret;
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        if (VERBOSE || verbose) {fprintf(stderr, "[debug::%s] > iter %d, treated %d\n", __func__, iter, ret);}
    }

    // clean up
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){
        stacku32_destroy(&buf[i]);
    }
    free(buf);
    stacku32_destroy(&candidates);
    if (VERBOSE || verbose) {
        fprintf(stderr, "[M::%s] treated %d\n", __func__, ret);
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

    for (v1=0; v1<auxsg->n_seq*2; v1++){
        if (hamt_asgarc_util_countSuc(auxsg, v1, 0, 0, base_label)!=2 || 
            hamt_asgarc_util_countPre(auxsg, v1, 0, 0, base_label)!=2
            ){
                continue;
            }
        if (hamt_asgarc_util_get_the_two_targets(auxsg, v1, &vtmp[0], &vtmp[1], 0, 0, base_label)!=0)
            continue;
        if (vtmp[0]!=v1 && vtmp[1]!=v1){continue;}  // v1 is not a self circle
        if (vtmp[0]==v1 and vtmp[1]==v1){  // sancheck, should not happen
            fprintf(stderr, "[E::%s] double arc at %.6d ? continue anyway\n", __func__, (int)(v1>>1));
            continue;
        }
        v0 = vtmp[0]==v1? vtmp[1] : vtmp[0];

        if (hamt_asgarc_util_get_the_two_targets(auxsg, v1^1, &vtmp[0], &vtmp[1], 0, 0, base_label)!=0)
            continue;
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

int hamt_ug_prectg_rescueShortCircuit(asg_t *sg, 
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
    int nb_modified = 0, nb_modified_tot=0;
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
        nb_modified_tot += nb_modified;
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
    return nb_modified_tot;
}

int hamt_ug_prectg_rescueShortCircuit_simpleAggressive(asg_t *sg, ma_ug_t *ug, 
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
    return nb_treat;
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
                       uint32_t *buf_blacklist,
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

// FUNC
//     Check all pairs of long unitigs/contigs.
//     If they can form a pair of source and sink,
//      try to find a path through.
// NOTE
//     Might be very slow given certain toplogy... 
//     Also ASan can be 7x slower.
int hamt_ug_resolveTangles(asg_t *sg, ma_ug_t *ug, 
                           int base_label, int alt_label){
    double startTime = Get_T();
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
            if (ug->u.a[vu>>1].subg_label==idx_subgraph){
                if (ug->u.a[vu>>1].len>70000) {has_long_tig = 1;}  // expecting at least one handle
                subgraph_size++;
                vecu32_push(&buf, vu);
                if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, base_label)!=0 &&
                    hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, base_label)!=0){
                        non_tip_tigs++;
                }
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
                    if (verbose) {fprintf(stderr, "[debug::%s] skip because either handle was a part of some tangles treated\n", __func__);}
                    continue;
                }

                if (verbose){
                    fprintf(stderr, "[debug::%s] check source %.6d sink %.6d \n", __func__, 
                                            (int)(source>>1)+1, (int)(sink>>1)+1);
                }

                // (note: must check label, since hamt_ug_arc_flip_between is used and it only flip labels for simplicity)
                if (base_label>=0 && auxsg->seq_vis[source>>1]!=base_label) {continue;}
                if (base_label>=0 && auxsg->seq_vis[sink>>1]!=base_label) {continue;}

                if (source==sink) {  // for hinted circular, it's source==(sink^1)
                    if (verbose) fprintf(stderr, "[deubg::%s] skip because sink==source\n", __func__);
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

    fprintf(stderr, "[M::%s] treated %d spots, used %.1fs\n", __func__, nb_treated, Get_T()-startTime);
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
        if (hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label)!=0)
            continue;
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
        if (hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label)!=0)
            continue;
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
    for (vu=0; vu<auxsg->n_seq*2; vu++){
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
            if (hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, base_label)!=0)
                continue;

            // get reads at the ends of wu0 and wu1
            if (wu[0]&1){
                startread = ug->u.a[wu[0]>>1].end^1;
            }else{
                startread = ug->u.a[wu[0]>>1].start;
            }
            if (wu[1]&1){
                endread = ug->u.a[wu[1]>>1].end^1;
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

    fprintf(stderr, "[M::%s] treated %d\n", __func__, cnt);
    if (cnt){
        free(sg->idx);
        sg->idx = 0;
        sg->is_srt = 0;
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }
    return cnt;
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
    //     Drop the shorter overlap at a bifurcation if the length diff between the two
    //       is extreme. Additionally, will consider:
    //       - coverage diff
    //       - topo: query has at least 2 non-tip targets (won't check further though)
    // RETURN
    //     number of overlaps dropped.
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
        if (hamt_asgarc_util_get_the_two_targets(auxsg, vu, &wu[0], &wu[1], 0, 0, 0)!=0)
            continue;
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
int hamt_ug_popLooseTangles_v2(asg_t *sg, ma_ug_t *ug, int max_step){  //TODO
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
    int verbose = 0;

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

void hamt_dump_haplotig_pairs(ma_ug_t *ug, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              char *asm_prefix){
    // FUNC
    //     Experimental; bin splitting assistance.
    //     Wrapper around "hamt_check_diploid_report2" (despite the name there's no diploid assumption).
    //     Write to tsv file: tigname1 tigname2 $hap_ratio
    char *output_name = (char*)malloc(strlen(asm_prefix)+30);
    char *utg_name = (char*)malloc(50);
    char *utg_name2 = (char*)malloc(50);
    sprintf(output_name, "%s.haplotigpairs.tsv", asm_prefix);
    FILE *fp = fopen(output_name, "w");
    assert(fp);

    asg_t *auxsg = ug->g;
    uint32_t vu, wu;
    float ratio[4];
    ma_utg_t *p;
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) continue;  // only need to check one direction
        for (wu=0; wu<auxsg->n_seq*2; wu++){  // yes check all pairs.
            if (wu&1) continue;
            p = &ug->u.a[vu>>1];
            sprintf(utg_name, "tig%.6d%c", (vu>>1) + 1, "lc"[p->circ]);
            p = &ug->u.a[wu>>1];
            sprintf(utg_name2, "tig%.6d%c", (wu>>1) + 1, "lc"[p->circ]);

            hamt_check_diploid_report2(ug, vu, wu, sources, &ratio[0], &ratio[1]);
            hamt_check_diploid_report2(ug, vu, wu, reverse_sources, &ratio[2], &ratio[3]);
            if (ratio[0]<0.0001 && ratio[1]<0.0001 && ratio[2]<0.0001 && ratio[3]<0.0001){
                continue;  // all (very close to) zero, no need to print.
            }
            fprintf(fp, "%s\t%d\t%s\t%d\t\%.4f\t%.4f\t%.4f\t%.4f\n", 
                    utg_name, (int)ug->u.a[vu>>1].len, utg_name2, (int)ug->u.a[wu>>1].len, 
                    ratio[0], ratio[1], ratio[2], ratio[3]);
            
        }
    }

    fclose(fp);
    free(output_name);
    free(utg_name); free(utg_name2);
}

int hamt_dump_path_coverage_with_haplotype_info_core(ma_ug_t *ug, asg_t *read_g, uint32_t vu, uint32_t wu, 
                                                      ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, 
                                                      R_to_U* ruIndex,
                                                      const ma_sub_t* coverage_cut,
                                                      uint64_t *buf, int buf_len, 
                                                      uint8_t *flag){
    // NOTE
    //     vu and wu both have the direction bit.
    //     `buf`: caller is responsible for allocating buf to length of at least utg_len(vu)+utg_len(wu)-overlapLength
    //     `flag`: caller is responsible to allocate to length of # reads. 
    //     (for the regular coverage est, see `get_ug_coverage`)
    // NOTE ON BIAS
    //     The estimations for ends, i.e. the start of vu and the end of wu, will be worse 
    //      than the middle segments. 
    // FUNC
    //     Given two connected unitigs/contigs vu and wu, 
    //      collected base-level infomation of cis- and trans- overlaps.
    //     Takes containments into consideration.
    //     For rescued/ditched/altered overlaps (e.g. hamt's containment treatment),
    //      their information is not retained, therefore this function will not account for that.
    // STORE RESULT
    //     `buf` is a packed counter. Each slot is a base. The packing is:
    //       cnt_het_total | cnt_hom_total | cnt_het | cnt_hom
    // RETURN
    //     0    fatal error (e.g. vu-wu arc not found)
    //     1    ok
    int verbose = 0;
    int san_not_el=0;

    asg_t *auxsg = ug->g;
    asg_arc_t ug_arc;
    if (!get_arc(auxsg, vu, wu, &ug_arc)){
        return 0;
    }
    // int buf_len = ug->u.a[vu>>1].len+ug->u.a[wu>>1].len-ug_arc.ol;
    assert (flag);
    assert (buf);
    memset(buf, 0, sizeof(uint64_t)*(buf_len) );
    memset(flag, 0, read_g->n_seq);

    ma_hit_t_alloc *the_source;
    int shift, offset_read, offset_tig;
    uint32_t is_unitig;
    uint32_t rId, rId_tn;
    uint8_t rev;
    ma_hit_t *h;
    ma_utg_t *utg_vu = &ug->u.a[vu>>1];
    if (ug->u.a[vu>>1].n==0 || ug->u.a[wu>>1].n==0) return 0;  // sancheck empty unitig (shouldn't happen)
    
    // count within the unitigs
    for (int which_utg=0; which_utg<=1; which_utg++){
        if (which_utg==0) {
            utg_vu = &ug->u.a[vu>>1];
            offset_tig = 0;
        }else {
            utg_vu = &ug->u.a[wu>>1];
            offset_tig = ug->u.a[vu>>1].len - ug_arc.ol;
            fprintf(stderr, "offset_tig is %d, vu len %d, ol %d\n", (int)offset_tig, (int)(ug->u.a[vu>>1].len), (int)(ug_arc.ol));
        }

        for (int which=0; which<=1; which++){
            if (which==0) {
                the_source = sources;
                shift = 0;
            }
            else {
                the_source = reverse_sources;
                shift = 16;
            }
            
            offset_read = 0;  // start position of read
            for (int k=0; k<utg_vu->n; k++){
                rId = utg_vu->a[k]>>33;
                
                // base covered by this read
                if (which==0){ 
                    flag[rId] = 1;
                    for (int i=0; i<coverage_cut[rId].e-coverage_cut[rId].s; i++){
                        // fprintf(stderr, "try to use %d\n", (int)(i+offset_tig+offset_read));
                        buf[i+offset_tig+offset_read] += ((uint64_t)2);  // 2 instead of 1 because i want a containment to count only as a half.
                    }
                }
                
                // also want to check if the current read contained any reads (that are now removed)
                // (similar to get_ug_coverage but without the r_flag)
                for (int i=0; i<the_source[rId].length; i++){
                    h = &the_source[rId].buffer[i];
                    if (h->el!=1) continue;
                    rId_tn = Get_tn((*h));
                    if (flag[rId_tn]) continue;
                    flag[rId_tn] = 1;
                    if (read_g->seq[rId_tn].del==1){
                        get_R_to_U(ruIndex, rId_tn, &rId_tn, &is_unitig);
                        if(rId_tn == (uint32_t)-1 || is_unitig == 1 || read_g->seq[rId_tn].del == 1) continue;
                    }
                    for (int i=Get_qs((*h)); i<Get_qe((*h)); i++){
                        buf[i+offset_tig+offset_read] += ((uint64_t)1)<<shift;  // containment counted as a half.
                    }
                }
                offset_read+=(uint32_t)utg_vu->a[k];  
                // fprintf(stderr, "(1)offset_read is %d\n", (int)offset_read);
            }
        }
        // end of: count within the unitigs

        // fprintf(stderr, "-------------------\n");

        // count outside of the unitigs
        for (int which_utg=0; which_utg<=1; which_utg++){
            if (which_utg==0) {
                utg_vu = &ug->u.a[vu>>1];
                offset_tig = 0;
            }else {
                utg_vu = &ug->u.a[wu>>1];
                offset_tig = ug->u.a[vu>>1].len - ug_arc.ol;
            }

            for (int which=0; which<=1; which++){
                if (which==0) {
                    the_source = sources;
                    shift = 0+32;
                }
                else {
                    the_source = reverse_sources;
                    shift = 16+32;
                }
                
                offset_read = 0;  // start position of read
                for (int k=0; k<utg_vu->n; k++){
                    rev = (utg_vu->a[k]>>32) & 1;
                    rId = utg_vu->a[k]>>33;
                    for (int i=0; i<the_source[rId].length; i++){
                        h = &the_source[rId].buffer[i];
                        /*if (h->el!=1) continue;*/  // commented out to allow inexact
                        rId_tn = Get_tn((*h));

                        if (flag[rId_tn]) continue;  // target read is in the contig, or contained by reads of the contig, or has already been counted in this section.

                        if (read_g->seq[rId_tn].del==1){
                            get_R_to_U(ruIndex, rId_tn, &rId_tn, &is_unitig);
                            if(rId_tn == (uint32_t)-1 || /*is_unitig == 1 ||*/ read_g->seq[rId_tn].del == 1) continue;  // allow the target read to be in other unitigs.
                        }
                        flag[rId_tn] = 1;
                        
                        // increment the counter: exact overlap
                        for (uint32_t j=(uint32_t)h->qns; j<h->qe; j++){
                            if (rev==0){
                                buf[j+offset_tig+offset_read] += ((uint64_t)2)<<shift;
                            }else{
                                buf[offset_tig+offset_read+read_g->seq[rId].len-j] += ((uint64_t)2)<<shift;
                            }
                        }
                        
                        // if there's any overhang, grab them into the counter of het with half the score
                        if ( ((uint32_t)h->qns==0 || h->qe>=read_g->seq[rId].len-1) && 
                             (h->ts==0 || h->te>=read_g->seq[rId_tn].len-1) ) continue;
                        for (int j=0; j<read_g->seq[rId_tn].len; j++){
                            if (j>=h->ts && j<h->te) continue;  // is the matched part
                            int tmp_offset;  // alignment block start position on the unitig (the position of qs if read is +, qe otherwise)
                            int tmp_loc;  // base position
                            if (rev==0){
                                tmp_offset = offset_tig + offset_read + (uint32_t)h->qns;
                            }else{
                                tmp_offset = offset_tig + offset_read + (read_g->seq[rId].len-h->qe);
                            }
                            if (h->rev==0){
                                tmp_loc = tmp_offset + j - h->ts;
                            }else{
                                tmp_loc = tmp_offset + h->te - j;
                            }
                            tmp_loc = tmp_loc>=0?  tmp_loc : 0;
                            tmp_loc = tmp_loc<buf_len? tmp_loc : buf_len-1;
                            // fprintf(stderr, "tmp_loc is %d\n", (int)(tmp_loc));
                            buf[tmp_loc] += ((uint64_t)1)<<48;
                        }
                        san_not_el++;
                        
                    }
                    offset_read+=(uint32_t)utg_vu->a[k]; 
                    // fprintf(stderr, "(2)offset_read is %d\n", (int)offset_read);
                }
            }
        }
        // endof: count outside of the unitigs
    }

    if (verbose){
        int san = 0;
        for (int i=0; i<read_g->n_seq; i++){
            if (flag[i]) san++;
        }
        fprintf(stderr, "[debug::%s] tigID vu %.6d wu %.6d, total reads considered %d (not el: %d), actual read counts were %d and %d\n", 
                            __func__,
                            (int)(vu>>1)+1, (int)(wu>>1)+1,
                            san, san_not_el, 
                            ug->u.a[vu>>1].n, ug->u.a[wu>>1].n);
    }

    return 1; // all ok
}
void hamt_dump_path_coverage_with_haplotype_info(ma_ug_t *ug, asg_t *read_g, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              R_to_U* ruIndex, const ma_sub_t* coverage_cut,
                              char *asm_prefix){
    // attempt of coverage estimation
    // writes down base level counts.
    char *output_name = (char*)malloc(strlen(asm_prefix)+30);
    char *utg_name = (char*)malloc(50);
    char *utg_name2 = (char*)malloc(50);
    uint8_t *flag = (uint8_t*)malloc(read_g->n_seq);
    int buf_m = 1000000, buf_n=0;
    uint64_t *buf = (uint64_t*)malloc(buf_m*sizeof(uint64_t));

    sprintf(output_name, "%s.hapcov.tsv", asm_prefix);
    FILE *fp = fopen(output_name, "w");
    assert (fp);

    asg_t *auxsg = ug->g;
    asg_arc_t *av;
    int nv;
    uint32_t vu, wu;
    float ratio[4];
    ma_utg_t *p;
    int shift_masks[4] = {48, 32, 16, 0};

    for (vu=0; vu<auxsg->n_seq*2; vu++){
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        if (nv==0) continue;
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            wu = av[i].v;
            p = &ug->u.a[vu>>1];
            sprintf(utg_name, "tig%.6d%c", (vu>>1) + 1, "lc"[p->circ]);
            p = &ug->u.a[wu>>1];
            sprintf(utg_name2, "tig%.6d%c", (wu>>1) + 1, "lc"[p->circ]);
            // fprintf(stderr, "> vu %s wu %s\n", utg_name, utg_name2);

            asg_arc_t ug_arc;
            if (!get_arc(auxsg, vu, wu, &ug_arc)){
                fprintf(stderr, "[E::%s] can't find wu %s (vu is %s)\n", __func__, utg_name2, utg_name);
            }

            if (vu==wu){
                buf_n = ug->u.a[vu>>1].len * 2;
                continue;  // won't choose to cut anyway, ignore it
            }
            else
                buf_n = ug->u.a[vu>>1].len+ug->u.a[wu>>1].len-ug_arc.ol;

            if (buf_n >buf_m){
                fprintf(stderr, "[debug::%s] buf expands from %d to %d\n", __func__, buf_m, buf_n);
                buf = (uint64_t*)realloc(buf, sizeof(uint64_t)*buf_n);
                buf_m = buf_n;
            }
            if (!hamt_dump_path_coverage_with_haplotype_info_core(ug, read_g, vu, wu,
                                                    sources, reverse_sources, 
                                                    ruIndex, coverage_cut, buf, buf_n, flag)){
                fprintf(stderr, "[E::%s] dump failed, ug1 %s ug2 %s\n", __func__, utg_name, utg_name2);
                continue;
            }
            fprintf(fp, "%s\t%s", utg_name, utg_name2);
            for (int i_mask=0; i_mask<4; i_mask++){
                fprintf(fp, "\t");
                for (int j=0; j<buf_n; j++){
                    fprintf(fp, "%d,", (int) ((uint16_t) (buf[j] >> shift_masks[i_mask]) ) )  ;
                }
            }
            fprintf(fp, "\n");
        }
    }


    fclose(fp);
    free(output_name);
    free(utg_name); free(utg_name2);
    free(flag);
    free(buf);
}


void hamt_update_coverage(ma_ug_t *ug, asg_t *read_g, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              R_to_U* ruIndex, const ma_sub_t* coverage_cut,
                              char *asm_prefix){
    // FUNC
    //     An expensive way to estimate coverage and coverage variances: not only count how many reads
    //      are there in the contig, but also consider all their overlapping reads, including trans-overlaps,
    //      at base level. 
    //     Will also write a tsv file.
    // NOTE
    //     Doesn't align well enough with alignment-based estimation..
    
    char *utg_name = (char*)malloc(50);
    char *output_name = (char*)malloc(strlen(asm_prefix)+30);
    sprintf(output_name, "%s.hapcov_mean_and_variance.tsv", asm_prefix);
    FILE *fp = fopen(output_name, "w");
    assert (fp);
    fprintf(fp, "contigName\tcontigLen\ttotalAvgDepth\t%s\t%s-var\n", asm_prefix, asm_prefix);

    asg_t *auxsg = ug->g;
    asg_arc_t *av;
    uint32_t nv;
    ma_utg_t *utg;
    ma_hit_t *h;
    uint32_t rId, rId_tn, is_unitig, uId, uId_tn;
    uint8_t rev;  // if the read is on reverse strand of a unitig   
    
    
    // mask of: which read could be used to estimate contig depth
    //   - non-pctg-contig read and isolated reads can be double-counted between contigs (not within one contig)
    //   - not-isolated reads that belong to the pctg (main assembly graph) will be only counted for its corresponding contig
    uint8_t *ban = (uint8_t*)malloc(read_g->n_seq);
    // init mask: collect reads that belong to the main assembly graph
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) continue;  // only need one direction
        if (hamt_asgarc_util_countPre(auxsg, vu, 0, 0, 0)==0 && hamt_asgarc_util_countSuc(auxsg, vu, 0, 0, 0)==0){
            continue;  // contig is unconnected, allow its read to be counted in other contigs
        }
        utg = &ug->u.a[vu>>1];
        for (int i=0; i<utg->n; i++){
            ban[utg->a[i]>>33] = 2;
        }
    }

    // count how many unitigs could a read overlap to
    // If a free-floating read (:=does not belong to any connected contigs
    //  in the primary assembly graph) have overlaps with more than one contigs,
    //  we don't want it to be counted as weight 1 every time.
    // Here implements not the best way, but better than ^ this: the weight
    //  is 1/nb_target_contigs. 
    uint8_t *targets = (uint8_t*)malloc(read_g->n_seq);
    for (uint32_t i=0; i<read_g->n_seq; i++){
        targets[i] = 1;
    }
    stacku32_t ru;
    stacku32_init(&ru);
    for (uint32_t i=0; i<sources->length; i++){
        if (ban[i]==2) continue;  // read belong to a contig, won't be counted more than one time anyway
        stacku32_reset(&ru);
        rId = i;
        get_R_to_U(ruIndex, rId, &uId, &is_unitig);
        for (uint32_t j=0; j<sources[i].length; j++){
            h = &sources[i].buffer[j];
            rId_tn = Get_tn((*h));
            get_R_to_U(ruIndex, rId_tn, &uId_tn, &is_unitig);
            if(rId_tn == (uint32_t)-1 || is_unitig == 1 || read_g->seq[rId_tn].del == 1) continue;
            if (uId_tn!=uId && !stack32_is_in_stack(&ru, uId_tn)){
                targets[rId] = targets[rId]==255? 255 : targets[rId]+1;
                stacku32_push(&ru, uId_tn);
            }
        }
    }

    // counter of per base depth
    int counter_n = 0;
    int counter_m = 5000000;
    float *counter = (float*)calloc(counter_m, sizeof(float));
    
    float mean, variance;
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (vu&1) continue;  // only need one direction
        if (ug->u.a[vu>>1].len>counter_m){
            counter_m = ug->u.a[vu>>1].len;
            counter = (float*)realloc(counter, sizeof(float)*counter_m);
        }
        counter_n = ug->u.a[vu>>1].len;
        utg = &ug->u.a[vu>>1];

        // reset mask and counter
        for (int i=0; i<read_g->n_seq; i++){
            if (ban[i]!=2) ban[i] = 0;
        }
        memset(counter, 0, sizeof(float)*counter_m);

        // collect counts
        int read_offset = 0;
        for (int k=0; k<utg->n; k++){
            rId = utg->a[k]>>33;
            // count the current read
            for (int idx=0; idx<coverage_cut[rId].e-coverage_cut[rId].s; idx++){
                counter[read_offset+idx]+=1;
            }

            // count any read that was contained by the current read
            for (int i=0; i<sources[rId].length; i++){
                h = &sources[rId].buffer[i];
                if (h->el!=1) continue;
                rId_tn = Get_tn((*h));
                if (ban[rId_tn]) continue;;
                if (read_g->seq[rId_tn].del==1){
                    get_R_to_U(ruIndex, rId_tn, &uId_tn, &is_unitig);
                    if(uId_tn == (uint32_t)-1 || is_unitig == 1 || read_g->seq[uId_tn].del == 1) continue;
                }
                for (int j=Get_qs((*h)); j<Get_qe((*h)); j++){
                    counter[read_offset+j]+= ((float)1)/targets[rId_tn];  // weighted
                }
                ban[rId_tn] = 1;
            }

            // count its overlapping reads
            rev = (utg->a[k]>>32)&1;
            for (int i=0; i<sources[rId].length; i++){
                h = &sources[rId].buffer[i];
                // if (h->el!=1) continue;  // do not allow inexact overlaps; comment out to allow.
                rId_tn = Get_tn((*h));
                if (ban[rId_tn]) continue;

                for (uint32_t j=(uint32_t)h->qns; j<h->qe;j++){
                    if (rev==0){
                        counter[read_offset+j]+=((float)1)/targets[rId_tn];
                    }else{
                        counter[read_offset+read_g->seq[rId].len-j]+=((float)1)/targets[rId_tn];
                    }
                }
                ban[rId_tn] = 1;
            }

            // update offset
            read_offset += (uint32_t)utg->a[k];
        }

        // calculate mean
        mean = 0;
        for (int i=0; i<counter_n; i++){
            mean += (counter[i]-mean)/(i+1);
        }

        // calculate std
        variance = 0;
        float tmp;
        for (int i=0; i<counter_n; i++){
            tmp = counter[i]-mean;
            variance += tmp*tmp;
        }
        variance = variance/counter_n;

        sprintf(utg_name, "s%d.ctg%.6d%c", (int)utg->subg_label, (vu>>1) + 1, "lc"[utg->circ]);
        fprintf(fp, "%s\t%d\t%.3f\t%.3f\t%.3f\n", utg_name, ug->u.a[vu>>1].len,
                                                    mean, mean, variance);
    }


    fclose(fp);
    free(output_name);
    free(utg_name);
    free(ban);
    free(counter);
    free(targets); stacku32_destroy(&ru);
}


/**
 * @brief Given the auxilliary graph of a bidirected graph, do SCC finding 
 *         and report the labels on the (ordinary) graph.\n 
 *        To get SCC labels of the bidirected grpah, use `hamt_ug_scc_test`.
 * @cite Ando, Kazutoshi, Satoru Fujishige, and Toshio Nemoto. 
 *        "Decomposition of a bidirected graph into strongly connected components and 
 *        its signed poset structure." 
 *        Discrete Applied Mathematics 68.3 (1996): 237-248.
 * @par g auxilliary graph of a bidirected graph
 * @par n_scc if available, will store # of SCCs
 * @return array of vertex SCC labels on the auxilliary graph, length is 
 * n_seq*2, must not contain -1.
*/
int* hamt_asggraph_scc(asg_t *g, int *n_scc){
    int verbose = 0;
    int *ret = (int*)malloc(sizeof(int)*g->n_seq*2);
    for (int i=0; i<g->n_seq*2; i++){
        ret[i] = -1;  // for sancheck
    }
    asg_t *gT = hamt_asggraph_util_gen_transpose(g);

    stacku32_t s;  // DFS stack
    stacku32_init(&s);
    uint8_t *color = (uint8_t*)calloc(g->n_seq*2, 1);
    uint64_t *end_time = (uint64_t*)calloc(g->n_seq*2, sizeof(uint64_t));  // packing: (time | vu)
    uint64_t time = 0;

    // collect timestamps
    uint32_t vu, wu;
    uint32_t nw;
    asg_arc_t *aw;
    int all_finished;
    for (vu=0; vu<g->n_seq*2; vu++){
        if (color[vu]) {
            if (color[vu]!=2){
                fprintf(stderr, "[E::%s] DFS bug (1)\n", __func__);
                exit(1);
            }
            continue;
        }
        // (DFS)
        stacku32_push(&s, vu);
        time++;
        color[vu] = 1;
        if (verbose) fprintf(stderr, "[debug::%s-1] seed utg%.6d\n", __func__, (int)(vu>>1)+1);
        while (stacku32_pop(&s, &wu)){
            if (verbose) fprintf(stderr, "[debug::%s-1]   pop utg%.6d\n", __func__, (int)(wu>>1)+1);
            if (color[wu]==2) continue;
            all_finished = 1;
            // check all child nodes
            nw = asg_arc_n(g, wu);
            aw = asg_arc_a(g, wu);
            for (int i=0; i<nw; i++){
                if (aw[i].del) continue;
                if (color[aw[i].v]==0){
                    if (all_finished){
                        if (verbose) fprintf(stderr, "[debug::%s-1]   put back utg%.6d\n", __func__, (int)(wu>>1)+1);
                        stacku32_push(&s, wu);
                        all_finished = 0;
                    }
                    stacku32_push(&s, aw[i].v);
                    if (verbose) fprintf(stderr, "[debug::%s-1]   push utg%.6d\n", __func__, (int)(aw[i].v>>1)+1);
                    color[aw[i].v] = 1;
                    time++;
                }
            }
            if (all_finished){
                if (verbose) fprintf(stderr, "[debug::%s-1]   finished utg%.6d\n", __func__, (int)(wu>>1)+1);
                color[wu] = 2;
                end_time[wu] = (time<<32) | ((uint64_t)wu);
                time++;
            }
            time++;
        }
    }

    radix_sort_ovhamt64(end_time, end_time+g->n_seq*2);

    // DFS on transposed graph
    memset(color, 0, g->n_seq*2);
    int scci = 0;
    for (int idx=g->n_seq*2-1; idx>=0; idx--){
        vu = (uint32_t)end_time[idx];
        if (color[vu]) {
            if (color[vu]!=2){
                fprintf(stderr, "[E::%s] DFS bug (2)\n", __func__);
                exit(1);
            }
            continue;
        }
        // (DFS)
        stacku32_push(&s, vu);
        color[vu] = 1;
        if (verbose) fprintf(stderr, "[debug::%s-2] seed utg%.6d\n", __func__, (int)(vu>>1)+1);
        while (stacku32_pop(&s, &wu)){
            if (verbose) fprintf(stderr, "[debug::%s-2]   pop utg%.6d\n", __func__, (int)(wu>>1)+1);
            if (color[wu]==2) continue;
            all_finished = 1;
            // check all child nodes
            nw = asg_arc_n(gT, wu);
            aw = asg_arc_a(gT, wu);
            for (int i=0; i<nw; i++){
                // (note: do not test aw[i].del here. It's no use && the transpose'd graph has no info.
                if (verbose) fprintf(stderr, "[debug::%s-2]      seeing utg%.6d\n", __func__, (int)(aw[i].v>>1)+1);
                if (color[aw[i].v]==0){
                    if (all_finished){
                        if (verbose) fprintf(stderr, "[debug::%s-2]   putback utg%.6d\n", __func__, (int)(wu>>1)+1);
                        stacku32_push(&s, wu);
                        all_finished = 0;
                    }
                    stacku32_push(&s, aw[i].v);
                    color[aw[i].v] = 1;
                    if (verbose) fprintf(stderr, "[debug::%s-2]   push utg%.6d\n", __func__, (int)(aw[i].v>>1)+1);
                }
            }
            if (all_finished){
                if (verbose) fprintf(stderr, "[debug::%s-2]   finished utg%.6d\n", __func__, (int)(wu>>1)+1);
                color[wu] = 2;
                ret[wu] = scci;
            }
        }
        scci++;
    }

    asg_destroy(gT);
    stacku32_destroy(&s);
    free(end_time);
    free(color);

    for (int i=0; i<g->n_seq*2; i++){
        if (ret[i]<0){
            fprintf(stderr, "[E::%s] some vertex not labelled\n", __func__);
            exit(1);
        }
    }

    if (n_scc) *n_scc = scci;
    return ret;
}

/**
 * @brief Label checker for bidirected graph SCC.
 * @par strict should be 1; 0 bypasses the XOR...
 * @return 1 if yes, 0 if no
*/
int hamt_ug_scc_test(int *labels, uint32_t v, uint32_t w, int strict){
    // FUNC
    //    Given SCC info from the aux graph of a bidirected graph, 
    //     test if the two vertices are strongly connected.
    //     (to generate the SCC status: hamt_asggraph_scc(auxsg, &n_scc0) )
    //    Ref: Ando et al 1996
    // RET
    //    0 if no
    //    1 if yes

    int label0, stat1, stat2;

    if (strict){
        label0 = labels[v<<1];
        stat1 = labels[(w<<1)| 0]==label0;
        stat2 = labels[(w<<1)| 1]==label0;
        if ( (stat1 && (!stat2)) || (stat2 && (!stat1)) ){
            return 1;
        }
        return 0;
    }else{
        return (labels[v<<1]==labels[w<<1] || 
                labels[v<<1]==labels[(w<<1) |1] ||
                labels[(v<<1) |1]==labels[w<<1] ||
                labels[(v<<1) |1]==labels[(w<<1) |1]);
    }
}
int hamt_ug_scc_count(ma_ug_t *ug, int *labels, uint32_t vu){
    // FUNC
    //     Count how many other vertices are strongly connected with vu (has direction)
    int ret = 0;
    asg_t *auxsg = ug->g;
    for (uint32_t wu=0; wu<auxsg->n_seq*2; wu++){
        if ((wu>>1)==(vu>>1)) continue;
        if (hamt_ug_scc_test(labels, vu>>1, wu>>1, 0)){
            ret++;
        }
    }
    return ret;
}

/**
 * @brief Given tig graph and SCC labels collected from aux graph, 
 * print tig and labels of the bidirected graph.
 * (Use hamt_utg_scc_debugprintext if you don't already have labels.)
 * @par ug unitig/contig graph
 * @par labels collected by `hamt_asggraph_scc`. Array length is 
 * n_seq*2.
*/
void hamt_utg_scc_debugprint(ma_ug_t *ug, int *labels){ 
    asg_t *auxsg = ug->g;
    int i_scc = 0;
    uint8_t *color = (uint8_t*)calloc(auxsg->n_seq, 1);
    for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
        if (color[vu>>1]) continue;
        fprintf(stderr, "[debug::%s]\tutg%.6d\t%d\n", __func__, 
                (int)(vu>>1)+1, i_scc);
        for (uint32_t wu=vu+2; wu<auxsg->n_seq*2; wu++){
            if (color[wu>>1]) continue;
            if ((labels[vu]==labels[wu]) ^ (labels[vu]==labels[wu^1])){
                fprintf(stderr, "[debug::%s]\tutg%.6d\t%d\n", __func__, 
                        (int)(wu>>1)+1, i_scc);
                color[wu>>1] = 1;
            }
        }
        color[vu>>1] = 1;
        i_scc++;
    }
}

/**
 * @brief For sprinkling into the main workflow. Give the tig graph, 
 * collect SCC labels and write to stderr. cut and awk friendly.
 * @par ug unitig/contig grpah.
*/
void hamt_utg_scc_debugprintext(ma_ug_t *ug){ 
    int n_scc0;
    asg_t *auxsg = ug->g;
    int *labels = hamt_asggraph_scc(auxsg, &n_scc0);
    fprintf(stderr, "%d\n", n_scc0);
    hamt_utg_scc_debugprint(ug, labels);
    free(labels);
}

#if 0
// ref: Donald B. Johnson (1975) Finding all the elementary circuits of a directed graph
// Helper function.
void hamt_ug_get_all_elementary_circuits_unblock(uint32_t vu, uint8_t *blocked, stacku32_t *B){
    // FUNC
    //     Helper routine to handle unblocking. 
    blocked[vu] = 0;
    uint32_t wu;
    for (int i=0; i<B[vu].n; i++){
        wu = B[vu].a[i];
        if (blocked[wu]) hamt_ug_get_all_elementary_circuits_unblock(wu, blocked, B);
    }
    stacku32_reset(&B[vu]);
}

// ref: Donald B. Johnson (1975) Finding all the elementary circuits of a directed graph
// Enumerate _all_ elemetary circuits.
int hamt_ug_get_all_elementary_circuits_circuit(ma_ug_t *ug, int *sgscc_labels,
                                                uint32_t root, uint32_t vu, 
                                                stacku32_t *stack, uint8_t *blocked, stacku32_t *B,
                                                stacku32_t *for_dump){
    int is_found = 0;
    int verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t nv, wu;
    asg_arc_t *av;

    stacku32_push(stack, vu);
    blocked[vu] = 1;
    nv = asg_arc_n(auxsg, vu);
    av = asg_arc_a(auxsg, vu);
    if (verbose) {fprintf(stderr, "[debug::%s] pushed utg%.6d\n", __func__, (int)(vu>>1)+1);}
    for (int i=0; i<nv; i++){
        wu = av[i].v;
        if ((wu>>1)<(root>>1)) continue;  // root has to be the smallest vertex
        if (verbose) {fprintf(stderr, "[debug::%s]  check wu=utg%.6d\n", __func__, (int)(wu>>1)+1);}
        if (!hamt_ug_scc_test(sgscc_labels, vu>>1, wu>>1, 0)) {
            if (verbose) fprintf(stderr, "[debug::%s]    vu-wu not SCC, skip wu\n", __func__);
            continue;
        } 
        
        if (wu==root) {
            // found a circuit, record it
            if (verbose) {
                fprintf(stderr, "[debug::%s]    found a cycle, stack size is %d, result stack size is %d\n", __func__, stack->n, for_dump->n);
                fprintf(stderr, "[debug::%s]      current stack: ", __func__);
                for (int j=0; j<stack->n; j++){
                    fprintf(stderr, "utg%.6d, ", (int)(stack->a[j]>>1)+1);
                }
                fprintf(stderr, "\n");
                fprintf(stderr, "[debug::%s]    current blocks: ", __func__);
                for (int j=0; j<auxsg->n_seq*2; j++){
                    if (blocked[j]) fprintf(stderr, "utg%.6d(%d), ", (int)(j>>1)+1, (int)(j&1));
                }
                fprintf(stderr, "\n");
            }
            
            fprintf(stdout, "ec-cycle\t");
            for (int j=0; j<stack->n; j++){
                // stacku32_push(for_dump, stack->a[j]);
                fprintf(stdout, "utg%.6d\t", (int)(stack->a[j]>>1)+1);
            }
            fprintf(stdout, "\n");
            for_dump->a[0]++;
            if (for_dump->a[0]>0xffffffff-2){
                fprintf(stderr, "[%s] too many circles\n", __func__);
                exit(1);
            }
            // stacku32_push(for_dump, (uint32_t)0xffffffff);  // delimiter
            is_found = 1;
        }else if (!blocked[wu]){
            if (verbose) fprintf(stderr, "[debug::%s]   wu not blocked, recursive call..\n", __func__);
            if (hamt_ug_get_all_elementary_circuits_circuit(ug, sgscc_labels, root, wu, stack, blocked, B, for_dump)){
                if (verbose) fprintf(stderr, "[debug::%s] returning from recurse, cycle found, current wu is utg%.6d\n", __func__, (int)(wu>>1)+1);
                is_found = 1;
            }
        }else{
            if (verbose) fprintf(stderr, "[debug::%s]   wu is blocked\n", __func__);
        }
    }

    if (is_found){
        hamt_ug_get_all_elementary_circuits_unblock(vu, blocked, B);
    }else{
        for (int i=0; i<nv; i++){
            wu = av[i].v;
            if ((wu>>1)<(root>>1)) continue;  // root has to be the smallest vertex
            blocked[wu] = 1;  // me: this is required, right??
            if (verbose) fprintf(stderr, "[debug::%s]    blocking wu=%.6d\n", __func__, (int)(wu>>1)+1);
            if (!stack32_is_in_stack(&B[wu], vu)){
                stacku32_push(&B[wu], vu);
                if (verbose) fprintf(stderr, "[debug::%s]    blocking vu=%.6d on wu=%.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
            }
        }
    }
    uint32_t tmp;
    stacku32_pop(stack, &tmp);
    assert(tmp==vu);
    return is_found;
}

// Enumerate _all_ elementary circuits.
// SCC and use Johnson's method
void hamt_ug_get_all_elementary_circuits(ma_ug_t *ug){
    int verbose = 0;
    asg_t *auxsg = ug->g;
    int n_scc0;
    int *labels = hamt_asggraph_scc(auxsg, &n_scc0);
    if (verbose) hamt_utg_scc_debugprint(ug, labels);

    exit(0);

    stacku32_t results, ec_stack;
    stacku32_init(&results);
    stacku32_init(&ec_stack);

    uint8_t *blocked = (uint8_t*)calloc(auxsg->n_seq*2, 1);
    stacku32_t *B = (stacku32_t*)malloc(auxsg->n_seq*2 * sizeof(stacku32_t));
    for (int i=0; i<auxsg->n_seq*2; i++){
        stacku32_init(&B[i]);
    }

    for (uint32_t root=0; root<auxsg->n_seq*2; root++){
        if (root&1) continue;  // only need to check in one direction
        if (hamt_ug_scc_count(ug, labels, root)==0) {  // no neighbor in SCC
            if (verbose) fprintf(stderr, "[debug::%s] root, skip utg%.6d because empty SCC\n", __func__, (int)(root>>1)+1);
            continue;  
        }
        if (verbose){
            fprintf(stderr, "[debug::%s] root=utg%.6d\n", __func__, (int)(root>>1)+1);
        }

        // reset blocking status
        memset(blocked, 0, auxsg->n_seq*2);
        for (int i=0; i<auxsg->n_seq*2; i++){
            stacku32_reset(&B[i]);
        }

        // search
        stacku32_push(&results, (uint32_t)0);
        hamt_ug_get_all_elementary_circuits_circuit(ug, labels, root, root, &ec_stack, blocked, B, &results);
    }

    // // debug
    // fprintf(stderr, "\nec-cycle\t");
    // for (int i=0; i<results.n; i++){
    //     if (results.a[i]==0xffffffff) fprintf(stderr, "\nec-cycle\t");
    //     else{
    //         fprintf(stderr, "tig%.6d\t", (int)(results.a[i]>>1)+1);
    //     }
    // }

    // cleanup
    free(blocked);
    stacku32_destroy(&results);
    stacku32_destroy(&ec_stack);
    for (int i=0; i<auxsg->n_seq*2; i++){
        stacku32_destroy(&B[i]);
    }
    free(B);
}
#endif


// FUNC
//     Given an array of unitig/contig IDs and the graph, 
// RETURN
//     Returns the sequence (heap allocation), 
//     with its length stored in seq_l.
char *hamt_ug_get_path_sequence(ma_ug_t *ug, uint32_t *r, int l, int is_circ, int *seq_l){
    int verbose = 0;
    int seq_m = 3000000, seq_n=0;
    char *seq = (char*)malloc(seq_m); 
    int seqtmp_m = 1000000;
    char *seqtmp = (char*)malloc(seq_m); 
    if (!seq || !seqtmp){
        fprintf(stderr, "[E::%s] malloc failed. seq_m=%d\n", 
                __func__, seq_m);
        if (seq) free(seq);
        if (seqtmp) free(seqtmp);
        *seq_l=0;
        return 0;
    }
   
    asg_t *auxsg = ug->g;
    ma_utg_t *p;
    char *s;
    int i, j, j2;
    uint32_t vu, wu;
    asg_arc_t *arc=0;
    asg_arc_t *av;
    uint32_t nv;
    int arclen=0, arclen_last=0;

    if (verbose) fprintf(stderr, "[debug::%s] start, cnt=%d\n", __func__, l);
    for (i=0; i<l; i++){
        vu = r[i];
        wu = r[i==l-1? 0 : i+1];

        p = &ug->u.a[vu>>1];
        if (vu&1){  // reverse strand
            if (p->len>=seqtmp_m){
                seqtmp_m = p->len+1;
                seqtmp = (char*)realloc(seqtmp, seqtmp_m); 
                if (!seqtmp){
                    fprintf(stderr, "[E::%s] realloc failed, seqtmp_m is %d\n", __func__, (int)seqtmp_m);
                    if (seq) free(seq);
                    if (seqtmp) free(seqtmp);
                    *seq_l = 0;
                    return 0;
                }
            }
            for (j=p->len-1, j2=0; j>=0; j--){
                seqtmp[j2] = seqcmp(p->s[j]);
                j2++;
            }
            seqtmp[j2] = '\0';
            s = seqtmp;
            if (verbose) {fprintf(stderr, "[debug::%s]   - %d ctg %.6d, rev, arclen %d, seqlen %d, seq_n %d, seq_m%d\n", 
                                __func__, i, (int)(vu>>1)+1, arclen, p->len, seq_n, seq_m); fflush(stderr);}
        }else{
            s = p->s;
            if (verbose) {fprintf(stderr, "[debug::%s]   - %d ctg %.6d, ori, arclen %d, seqlen %d, seq_n %d, seq_m%d\n", 
                                __func__, i, (int)(vu>>1)+1, arclen, p->len, seq_n, seq_m); fflush(stderr);}
        }

        // push the contig (if it's not the last one)
        if (i!=l-1){
            if (seq_n + p->len - arclen>=seq_m){
                seq_m += p->len - arclen;
                seq = (char*)realloc(seq, seq_m); 
                if (!seq){
                    fprintf(stderr, "[E::%s] realloc failed-b, seq_m is %d\n", __func__, (int)seq_m);
                    if (seq) free(seq);
                    if (seqtmp) free(seqtmp);
                    *seq_l = 0;
                    return 0;
                }
            }
            sprintf(seq+seq_n, "%.*s", (int)(p->len-arclen), s+arclen);
            seq_n += p->len-arclen;
        }
        else
            arclen_last = arclen;
        
        // update arc length 
        av = asg_arc_a(auxsg, vu);
        nv = asg_arc_n(auxsg, vu);
        for (j=0; j<nv; j++){
            if (av[j].v==wu){
                arc = &av[j];
                break;
            }
        }
        if ((i==l-1 && is_circ) || i<l-1){  // must find the arc, unless it's the end of a non-circular path
            assert(j<nv);
            arclen = (int)arc->ol;
        }else{
            arclen = 0;
        }

        // push the last contig: it is truncated on both ends if path is circular,
        //                       or, the collected sequence needs a truncation.
        if (i==l-1){
            if (seq_n + p->len - arclen - arclen_last>=seq_m){
                seq_m += p->len - arclen - arclen_last;
                seq = (char*)realloc(seq, seq_m); 
                if (!seq){
                    fprintf(stderr, "[E::%s] realloc failed-c, seq_m is %d\n", __func__, (int)seq_m);
                    if (seq) free(seq);
                    if (seqtmp) free(seqtmp);
                    *seq_l = 0;
                    return 0;
                }
            }
            if (p->len-arclen-arclen_last>=0){
                if (verbose) {fprintf(stderr, "[debug::%s]   - last status: type add, len=%d, arclen=%d, arclen_last=%d\n", 
                                                __func__, p->len, arclen, arclen_last); fflush(stderr);}
                sprintf(seq+seq_n, "%.*s", (int)(p->len-arclen-arclen_last), s+arclen_last);
                seq_n += p->len - arclen - arclen_last;
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s]   - last status: type shrink, len=%d, arclen=%d, arclen_last=%d\n", 
                                                __func__, p->len, arclen, arclen_last); fflush(stderr);}
                seq_n -= arclen + arclen_last - p->len;
            }
        }
    }

    free(seqtmp);
    *seq_l = seq_n;
    return seq;
}

typedef struct {
    // input
    ma_ug_t *ug; 
    stacku32_t *paths;
    int *offsets;  
    // output
    char **seqs;
    int *seqs_ll;
}hamt_getseq_t; 

static void hamt_ug_get_path_sequence_parallel_callback(void *data, long cid, int tid){
    hamt_getseq_t *d = (hamt_getseq_t*)data;
    int offset = d->offsets[cid];
    stacku32_t *s = d->paths;
    d->seqs[cid] = hamt_ug_get_path_sequence(d->ug, s->a+offset+1, s->a[offset], 1, &d->seqs_ll[cid]);
}

/**
 * @func Paralleled hamt_ug_get_path_sequence, store results.
 * @par seqs_ll_p Stores outputs, seq lengths.
 * @par seqs_p Output, seqs, newly allocated.
*/
void hamt_ug_get_path_sequence_parallel(ma_ug_t *ug, stacku32_t *s, 
                                          int n_paths, int n_threads, 
                                          int *seqs_ll_p, char **seqs_p){
    // collect
    hamt_getseq_t *d = (hamt_getseq_t*)calloc(1, sizeof(hamt_getseq_t));
    d->ug = ug;
    d->paths = s;
    d->offsets = (int*)calloc(n_paths, sizeof(int)); assert(d->offsets);
    d->seqs = seqs_p;
    d->seqs_ll = seqs_ll_p;

    int offset=0, idx=0, l=0;
    int seq_l;   
    while (idx<n_paths){//while (offset<s->n){
        if (offset>=s->n){
            fprintf(stderr, "[E::%s] offset OOB reading, check the path stack. Aborting.\n", __func__);
            fprintf(stderr, "[E::%s] n_paths=%d, dump of path stack:\n", __func__, n_paths);
            for (int i=0; i<s->n; i++){
                fprintf(stderr, "[dump::%s] i=%d val=%d\n", __func__, i, (int)s->a[i]);
            }
            exit(1);
        }
        d->offsets[idx] = offset;
        l = s->a[offset];  // packed as: [n_tigs, ID_tig1, ID_tig2, ...., ID_tign, n_tigs2, ID_tig1, ...]
        offset += l+1;
        idx++;
    }

    kt_for(n_threads, hamt_ug_get_path_sequence_parallel_callback, d, n_paths);

    // cleanup
    free(d->offsets);
    free(d);
}

// opportunistic elementary circuits: 
//   - seed from long contigs (TODO)
//   - global soft block:  don't extend if one try uses too many "used" contigs
//   - DFS stepping favors "unused" white contig first
//   - unused contig contig gives a minus weight instead of zero,
//     to encourage collecting any path that (somewhat) supersets 
//     a shorter, established path.
int hamt_asg_get_one_cycle_with_constraint(ma_ug_t *ug, uint32_t root, 
                                                         int *scclables, double *weights, uint8_t *color0,
                                                         int min_length, 
                                                         int max_weight,
                                                         stacku32_t *report_stack){
    int ret=0;
    int verbose = 0;  // print DFS debug info
    int is_print_found = 0;  // print out the cycle
    asg_t *g = ug->g;
    uint32_t nv, vu, wu;
    asg_arc_t *av;
    stacku32_t s, sb;
    stacku64_t scoring;
    uint8_t *color;
    int exhausted;
    if (color0) {color = color0; }
    else color = (uint8_t*)calloc(g->n_seq*2, sizeof(uint8_t));

    // constraint counters
    int total_length = 0;
    int total_weight = 0;
    int n_drop, i_d;

    stacku32_init(&s);
    stacku32_init(&sb);
    stacku64_init(&scoring);
    stacku32_push(&s, root);
    total_length += ug->u.a[root>>1].len;
    total_weight += weights[root>>1]<1? -1 : weights[root>>1];
    color[root] = 1;
    while (stacku32_pop(&s, &vu)){
        if (verbose) {fprintf(stderr, "[debug::%s] pop %.6d\n", __func__, (int)(vu>>1)+1);}
        if (color[vu]==2) continue;

        nv = asg_arc_n(g, vu);
        av = asg_arc_a(g, vu);
        stacku64_reset(&scoring);
        for (int i=0; i<nv; i++){
            stacku64_push(&scoring, ((uint64_t)(weights[av[i].v>>1]*100))<<32 | (uint64_t)i );
        }
        radix_sort_ovhamt64(scoring.a, scoring.a+scoring.n);
        exhausted = 1;
        for (int i=0; i<nv; i++){
            wu = av[(uint32_t)scoring.a[i]].v;
            if (wu==root){
                if (total_length>=min_length){
                    ret = 1;
                    stacku32_push(&s, vu);

                    // // // one more try: truncate 25%, 50% and 75% and try again, can we find a smaller circle?
                    // if (verbose) {
                    //     fprintf(stderr, "[debug::%s] truncation try:::\n", __func__);
                    // }
                    // stacku32_reset(&sb);
                    // stacku32_copyover(&s, &sb);
                    // for (float r=0.25; r<=0.75; r+=0.25){
                    //     if (verbose) {fprintf(stderr, "[debug::%s] truncation try::: r=%f\n", __func__, r);}
                    //     n_drop=s.n/4;
                    //     if (n_drop==0) break;
                    //     for(int j=0; j<n_drop; j++){
                    //         stacku32_pop(&s, &wu);
                    //         color[wu] = 2;
                    //     }
                    //     hamt_asg_get_one_cycle_with_constraint(ug, root, scclables, weights, color0, 
                    //                                             min_length, max_weight, report_stack);
                    // }
                    // stacku32_reset(&s);
                    // stacku32_copyover(&sb, &s);

                    // push the circuit
                    if (is_print_found) fprintf(stderr, "[debug::%s] cycle(root=utg%.6d), unitigs' total length=%d:\n", __func__, (int)(root>>1)+1, total_length);
                    stacku32_push(report_stack, s.n);
                    for (int j=0; j<s.n; j++){
                        if (is_print_found) fprintf(stderr, "[debug::%s]         utg%.6d%c\n", __func__, (int)(s.a[j]>>1)+1, "+-"[s.a[j]&1]);
                        stacku32_push(report_stack, s.a[j]);
                    }
                    // update weights
                    for (int j=0; j<s.n; j++){
                        // weights[s.a[j]>>1] += log10( (double) (ug->u.a[s.a[j]>>1].len) );  // length-dependent weight
                        weights[s.a[j]>>1] ++;  // a simple weight
                    }
                    goto finish;
                }
            }
            if (!hamt_ug_scc_test(scclables, vu>>1, wu>>1, 0)) continue;
            if (color[wu]) continue;
            
            // extra constraints: check length
            if (/*total_length+ug->u.a[wu>>1].len > max_length ||*/ total_weight+weights[wu>>1] > max_weight) {
                // Instead of killing this search and continue DFS,
                //   try to drop the last 1/4 of the DFS stack (marking them as visited) to force a retry.
                // Will not mark child nodes of these dropped nodes as visited, though.
                if (verbose) {fprintf(stderr, "[debug::%s] ======weight constraint failed, chop and retry\n", __func__);}
                if (s.n<4) goto finish;
                n_drop = s.n/4;
                exhausted = 1;
                for (i_d=0; i_d<n_drop; i_d++){
                    stacku32_pop(&s, &wu);
                    color[wu] = 2;
                    total_weight -= weights[wu>>1];
                }
                if (verbose) {fprintf(stderr, "[debug::%s] ======search head now: %.6d\n", __func__, (int)(s.a[s.n-1]>>1)+1);}
                break;  // leave the push loop
            }

            // (the regular white node push)
            if (verbose) {fprintf(stderr, "[debug::%s]   put back %.6d, push %.6d\n", __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);}
            exhausted = 0;
            stacku32_push(&s, vu);
            stacku32_push(&s, wu);
            total_length += ug->u.a[wu>>1].len;  // TODO: should minus overlap length! But it's expensive to get. Maybe just arbitrarily minus 8k or something.
            total_weight += weights[root>>1]<1? -1 : weights[wu>>1];
            color[wu] = 1;
            break;
            
        }
        if (exhausted) {
            color[vu] = 2;
            total_length -= ug->u.a[vu>>1].len;  // TODO: same, shouldn't be the whole unitig length
            total_weight -= weights[vu>>1];
        }
    }


finish:
    stacku32_destroy(&s);
    stacku32_destroy(&sb);
    stacku64_destroy(&scoring);
    if (!color0) free(color);
    return ret;

}
int hamt_ug_sum_lengths_of_unitigs(uint32_t *a, int n, ma_ug_t *ug, int is_path, int is_circ){
    // FUNC
    //     Given an array of unitig IDs (has direction but will ignore), 
    //      return a rough estimation of path length if it's a path,
    //      sum of lengths otherwise.
    //     If input is a path, the array is expected to be in order.
    // PAR
    //     a - array of length n, containing unitig IDs (with the direction bit)
    int ret = 0;
    for (int i=0; i<n; i++){
        ret += ug->u.a[a[i]>>1].len;
    }
    if (is_path){  // minus overlap lengths
        uint32_t nv, vu, wu;
        asg_arc_t *av;
        asg_t *auxsg = ug->g;
        for (int j, i=0; i<n-1; i++){
            vu = a[i];
            wu = a[i+1];
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            for (j=0; j<nv; j++){
                if (av[j].v!=wu) continue;
                ret-=(uint32_t)av[j].ul;
                break;
            }
            if (j==nv){
                fprintf(stderr, "[E::%s] looking for overlap, not found (linear); continue anyway\n", __func__);
            }
        }
        if (is_circ){
            vu = a[n-1];
            wu = a[0];
            nv = asg_arc_n(auxsg, vu);
            av = asg_arc_a(auxsg, vu);
            int j;
            for (j=0; j<nv; j++){
                if (av[j].v!=wu) continue;
                ret-=(uint32_t)av[j].ul;
                break;
            }
            // if (j==nv){
            //     fprintf(stderr, "[E::%s] looking for overlap, not found (circ link %.6d -> %.6d); continue anyway\n", 
            //             __func__, (int)(vu>>1)+1, (int)(wu>>1)+1);
            // }
        }
    }
    return ret;
}


// FUNC
//     Given a stack consists of found paths, drop paths that are too similar to 
    //     Given a stack consists of found paths, drop paths that are too similar to 
//     Given a stack consists of found paths, drop paths that are too similar to 
//      others (by checking unitig ID or haplotig status).
//     WILL modify `s` in-place.
// NOTE
//     `s` is filled like so: array length1, value1, value2...., array length2, value1, value2, ... 
//     Long contig is always preferred over path-finding.  
void hamt_ug_opportunistic_elementary_circuits_helper_deduplicate(stacku32_t *s, ma_ug_t *ug){
    int verbose = 0;
    int of1=0, of2=0, of1idx=0, of2idx=1;
    int hits=0, hits_len=0;  // count duplications
    int ll1, ll2;  // total unitig/contig lengths in path1 and path2, for checking duplication ratio in base pairs
    uint32_t l1, l2;
    int stat;
    if (s->n==0){
        fprintf(stderr, "[W::%s] empty stack\n", __func__);
        return ;
    }
    of2 = s->a[0]+1;

    stacku32_t cmp;
    stacku32_init(&cmp);

    // get number of paths in the stack, and length of each path
    uint8_t *color;
    // uint64_t *pathlengths;
    {
        int l=0;  // number of paths
        while (of1<s->n){
            l++;
            of1+=s->a[of1]+1;
        }
        color = (uint8_t*)calloc(l, 1);  // makr whether a path has been deleted
        // pathlengths = (uint64_t*)calloc(l, sizeof(uint64_t));
        // of1 = 0;
        // int i=0; 
        // while (of1<s->n){
        //     stacku32_reset(&cmp);
        //     l1 = s->a[of1];
        //     for (int i=0; i<l1; i++){
        //         stacku32_push(&cmp, s->a[of1+1+i]);
        //     }
        //     pathlengths[i++] = ((uint64_t)hamt_ug_sum_lengths_of_unitigs(cmp.a, l1, ug, 1, 1))<<32 | (uint32_t)i;
        // }
        // radix_sort_ovhamt64(pathlengths, pathlengths+l);
    }


    
    uint32_t maxl1, maxl2;  // length of the longest contig of a path
    of1 = 0;
    while (of1<s->n){
        maxl1 = 0;
        maxl2 = 0;

        // collect path1
        l1 = s->a[of1];
        if (color[of1idx]){  // path masked by a previous search, no need to check for more matches.
            of1idx++;
            of1+=l1+1;
            continue;
        }
        stacku32_reset(&cmp);
        for (int i=0; i<l1; i++){
            stacku32_push(&cmp, s->a[of1+1+i]);
            maxl1 = maxl1>s->a[of1+1+i]? maxl1 : s->a[of1+1+i];
        }
        ll1 = hamt_ug_sum_lengths_of_unitigs(cmp.a, l1, ug, 1, 1);
        if (ll1>8000000 || ll1<1000000) {  // path too large or too small
            color[of1idx] = 1;
            if (verbose){
                fprintf(stderr, "[debug::%s] dropping1=%c (size;%d)\n", __func__, of1idx, ll1);
            }
            continue;
        }
        
        of2 = of1+1+l1;
        of2idx = of1idx+1;

        while (of2<s->n){
            l2 = s->a[of2];
            for (int i=0; i<l2; i++) {
                maxl2 = maxl2>s->a[of2+1+i]? maxl2 : s->a[of2+1+i];
            }
            
            // collect path2
            ll2 = hamt_ug_sum_lengths_of_unitigs(&s->a[of2+1], l2, ug, 1, 1);
            if (ll2>8000000 || ll2<1000000){  // path too large
                color[of2idx] = 1;
                if (verbose){
                    fprintf(stderr, "[debug::%s] dropping2=%c (size;%d)\n", __func__, of2idx, ll2);
                }
                continue;
            }
            // (do nothing if both paths have long contigs)
            if (maxl1>300000 && maxl2>300000 && (maxl1!=maxl2)) {
                if (verbose){
                    fprintf(stderr, "[debug::%s] at %c, did not check %d (size)\n", 
                        __func__, of1idx, of2idx);
                }
                continue;
            }
            // guess similarity: path 2 compared to path 1
            hits_len = hits = 0;
            for (int j=0; j<l2; j++){
                uint32_t wu = s->a[of2+1+j];
                stat = stacku32_index_value(&cmp, wu);
                if (stat>=0){
                    hits++;
                    hits_len += ug->u.a[wu>>1].len;
                }else{// otherwise, check if it's a haplotig pair
                    for (int i=0; i<l1; i++){
                        stat = hamt_check_diploid(ug, cmp.a[i], wu, 0.5, R_INF.reverse_paf);
                        if (stat>0){
                            hits++;
                            hits_len += ug->u.a[wu>>1].len;
                        }
                    }
                }
            }

            // do we want to mask path2? Or path1?
            if ( /*(float)hits/(float)(l1>l2? l2 : l1) > 0.5 ||*/
                 (float)hits_len/(float)(ll1>ll2? ll2 : ll1) > 0.5 ){
                if (!(maxl2>300000 && maxl1<300000)) {
                    color[of2idx] = 1;
                    // self note: don't leave the loop here to remove the most duplicates.
                    //            You can leave and (thus randomly) let some cycles survive, though. 
                }
            }
            if (verbose){
                fprintf(stderr, "[debug::%s] dropping=%c: ", __func__, "NY"[color[of2idx]]);
                for (int i=0; i<l2; i++){
                    fprintf(stderr, "%.6d, ", (int)(s->a[of2+1+i]>>1)+1);
                }
                fprintf(stderr, "\n[debug::%s]   base: ", __func__);
                for (int i=0 ;i<cmp.n; i++){
                    fprintf(stderr, "%.6d, ", (int)(cmp.a[i]>>1)+1);
                }
                fprintf(stderr, "\nhits=%d, hits_len=%d, ll1=%d, ll2=%d\n", hits, hits_len, ll1, ll2);
            }

            // step
            of2+=l2+1;
            of2idx++;
        }
        of1 += l1+1;
        of1idx++;
    }

    stacku32_reset(&cmp);
    for (of1=of1idx=0; of1<s->n; of1idx++){
        l1 = s->a[of1];
        if (!color[of1idx]){
            stacku32_push(&cmp, (uint32_t)l1);
            for (int i=0; i<l1; i++){
                stacku32_push(&cmp, s->a[of1+1+i]);
            }
        }
        of1+=l1+1;
    }
    // overwrite
    stacku32_reset(s);
    stacku32_copyover(&cmp, s);

    free(color);
    stacku32_destroy(&cmp);
}


// FUNC
//     Bottom minhash-based sequence divergence estimation and deduplication.
// RETURN
//     Number of remaining paths.    
int hamt_ug_opportunistic_elementary_circuits_helper_deduplicate_minhash(stacku32_t *s, ma_ug_t *ug, 
                                                                        int kmersize, int n_hash,
                                                                        int n_thread){
    if (kmersize>31){
        fprintf(stderr, "[W::%s] k (%d) too large, reduced to 31\n", __func__, kmersize);
        kmersize = 31;
    }
    if (kmersize<3){
        fprintf(stderr, "[W::%s] k (%d) too small, will use default instead.\n", __func__, kmersize);
        kmersize = 21;
    }

    int verbose = 0;
    int n_paths = 0, n_remains=0;
    double time = Get_T();
    uint8_t *mask;  // marks whether a path has been deleted
    
    // get the number of paths in the stack
    {
        int offset = 0;
        while (offset<s->n){
            n_paths++;
            offset += s->a[offset]+1;
        }
        mask = (uint8_t*)calloc(n_paths, 1); assert(mask);
    }
    if (verbose) fprintf(stderr, "[debug::%s] got %d paths\n", __func__, n_paths);

    // collect sequences
    time = Get_T();
    int *seqs_ll = (int*)calloc(n_paths, sizeof(int)); assert(seqs_ll);
    char **seqs = (char**)malloc(sizeof(char*)*n_paths); assert(seqs);
    hamt_ug_get_path_sequence_parallel(ug, s, n_paths, n_thread, seqs_ll, seqs);
    fprintf(stderr, "[T::%s] got the sequences, used %.1fs\n", __func__, Get_T()-time);

    // sketch and pariwise
    time = Get_T();
    double **dists = hamt_minhash_mashdist(seqs, seqs_ll, n_paths, kmersize, n_hash, n_thread);
    fprintf(stderr, "[T::%s] collected mash distances for %d seqs, used %.1fs\n", 
                    __func__, n_paths, Get_T()-time); fflush(stderr);
    if (verbose) {
        for (int i=0; i<n_paths; i++){
            for (int j=i+1; j<n_paths; j++){
                fprintf(stderr, "[debug::%s] %d\t%d\t%.6f\n", __func__, i, j, dists[i][j]);
            }
        }
    }


    // drop some
    time = Get_T();
    int n_skip_lengthdiff=0;
    int n_skip_lengthabs=0;
    for (int i=0; i<n_paths; i++){
        if (mask[i]) continue;
        if (seqs_ll[i]<1000000 || seqs_ll[i]>8000000){  // weird size, drop the candidate
            mask[i] = 1;
            n_skip_lengthdiff++;
            continue;
        }
        for (int j=i+1; j<n_paths; j++){
            if (mask[j]) continue;
            if (seqs_ll[j]<1000000 || seqs_ll[j]>8000000){
                mask[j] = 1;
                n_skip_lengthdiff++;
                continue;
            }
            // do nothing if large length difference between the two paths 
            if ((seqs_ll[i]>seqs_ll[j]? seqs_ll[i]-seqs_ll[j] : seqs_ll[j]-seqs_ll[i])>1000000 ){
                n_skip_lengthdiff++;
                continue;
            }
            // check mash distance
            if (dists[i][j]<0.05)
                mask[j] = 1;
        }
    }

    // collect the remaining paths
    stacku32_t remain;
    stacku32_init(&remain);
    int of1, idx1, l1;
    for (of1=0, idx1=0; of1<s->n; idx1++){
        l1 = s->a[of1];
        if (mask[idx1]){
            of1 += l1+1;
            continue;
        }
        for (int i=of1; i<of1+l1+1; i++){
            stacku32_push(&remain, s->a[i]);
        }
        of1 += l1+1;
        n_remains++;
    }
    stacku32_reset(s);
    stacku32_copyover(&remain, s);

    stacku32_destroy(&remain);
    for (int i=0; i<n_paths; i++){
        free(dists[i]);
        free(seqs[i]);
    }
    free(dists);
    free(seqs);
    free(seqs_ll);
    free(mask);

    fprintf(stderr, "[M::%s] had %d paths, %d remained (%d dropped by length diff, %d by length abs)," 
                    "used %.1fs after sketching.\n", __func__, 
                    n_paths, n_remains,
                    n_skip_lengthdiff, n_skip_lengthabs, Get_T()-time);
    return n_remains;
}



// FUNC
//     Given a list of unitig/contig IDs, write the path's sequence in fasta format.
//     Will not check whether the path is valid; caller is responsible.
// PAR
//     idx is used for naming the resulting contig.
//     s is the start index. e is the end index (non-inclusive).
void hamt_ug_write_path_to_fasta(ma_ug_t *ug, uint32_t *r, int l, int is_circ, 
                                int pathID, FILE *fp){
    int i, j, j2=0, i2;
    int seq_l;
    char *seq;

    uint32_t vu, wu;
    asg_arc_t *arc;
    asg_arc_t *av;
    uint32_t nv;
    int arclen=0, arclen_last;

    fprintf(fp, ">rc.%d", pathID);
    for (i=0; i<l; i++){
        fprintf(fp, " s%d.ctg%.6d%c%c", (int)ug->u.a[r[i]>>1].subg_label, 
                (int)(r[i]>>1)+1, "lc"[ug->u.a[r[i]>>1].circ],
                "-+"[r[i]&1]);
    }
    fprintf(fp, "\n");
    
    seq = hamt_ug_get_path_sequence(ug, r, l, is_circ, &seq_l);
    fprintf(fp, "%.*s\n", seq_l, seq);
    free(seq);
}


// FUNC
//     Report deduplicated elementary circuits from a given unitig graph.
//     This function does not try to enumerate all combinations
//      as there are too many; instead it tries to get as many different circuits
//      as possible, then deduplicate using read overlap information.
// RET
//     List of long contigs >=100kb used in the deduplicated circles.
//     To be consumed by binning functions.
// PRE-CONDITION NOTE
//     Unitig graph's sequence must be collected.   
vu32_t *hamt_ug_opportunistic_elementary_circuits(asg_t *sg, ma_ug_t *ug, int n_thread){
    double time = Get_T();
    int is_print_found = 0;  // 0 to print, 1 to write a p_ctg2 file.
    FILE *fp = 0;
    char *fname = (char*)malloc(strlen(asm_opt.output_file_name)+50);
    
    asg_t *auxsg = ug->g;
    int n_scc0;
    int *labels = hamt_asggraph_scc(auxsg, &n_scc0);
    // (debug print: output SCC label - note that testing if a pair of vertices are
    //  in the same SCC is not simply comparing their labels.)
    // fprintf(stderr, ">>>\n");
    // for (uint32_t vu=0; vu<auxsg->n_seq*2; vu++){
    //     // if (vu&1) continue;
    //     fprintf(stderr, "SCCdump\tctg%.6d%c\t%d\n", (int)(vu>>1)+1, "+-"[vu&1], labels[vu]);
    // }
    
    int total_report = 0, prev_total_report=-1;
    int ret = 0;
    double *weights = (double*)calloc(auxsg->n_seq, sizeof(double));
    uint8_t *color = (uint8_t*)malloc(auxsg->n_seq*2);
    stacku32_t report_stack;
    stacku32_init(&report_stack);

    // stuff used to skip some searches
    uint8_t *used = (uint8_t*)calloc(auxsg->n_seq, 1);  // if a vertex has been used in a cycle, never use it as root
    int idx=0;

    int *scc_counts = (int*)calloc(auxsg->n_seq, sizeof(int));
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){
        if (i&1) continue;
        scc_counts[i>>1] = hamt_ug_scc_count(ug, labels, i);
    }
    
    while (prev_total_report<total_report){
        prev_total_report = total_report;
        for (uint32_t root=0; root<auxsg->n_seq*2; root++){
            if (root&1) continue;
            if (scc_counts[root>>1]==0) continue;  // vertex is on its own.
            if (used[root>>1]) continue; // vertex is already in a cycle

            memset(color, 0, sizeof(uint8_t)*auxsg->n_seq*2);
            
            ret = hamt_asg_get_one_cycle_with_constraint(ug, root, labels, weights, color, 
                                                        500000, 100, &report_stack);
            if (ret) { // found a new cycle
                fprintf(stderr, "[dbg::%s] cyc from idx %d to %d in report stack\n", 
                        __func__, idx, report_stack.n);
                total_report+=1;
                // udpate mask
                int tmpi=0, tmpl;
                uint32_t *tmpa = report_stack.a+idx;
                while (tmpi<report_stack.n-idx){
                    tmpl = tmpa[tmpi];
                    for (int i=tmpi; i<tmpl; i++){
                        used[tmpa[i+1]>>1] = 1;
                    }
                    tmpi+=tmpl+1;
                }
                idx = report_stack.n;
            }
        }
    }
    fprintf(stderr, "[M::%s] collected %d circuits, used %.2fs\n", __func__, total_report, Get_T()-time);
    
    int of=0, l=0, ofidx=0;

    // dump all found circuits, no deduplication
    time = Get_T();
    if (!is_print_found){
        sprintf(fname, "%s.rescue.all.fa", asm_opt.output_file_name);
        fp = fopen(fname, "w"); assert(fp);
        while (of<report_stack.n){
            l = report_stack.a[of];
            hamt_ug_write_path_to_fasta(ug, report_stack.a+of+1, l, 1, ofidx, fp);
            of = of+1+l;
            ofidx++;
        }
        fclose(fp);
    }
    fprintf(stderr, "[M::%s] wrote all rescued circles, used %.2fs\n", __func__, Get_T()-time);

    // deduplication and dump
    time = Get_T();
    // hamt_ug_opportunistic_elementary_circuits_helper_deduplicate(&report_stack, ug);
    hamt_ug_opportunistic_elementary_circuits_helper_deduplicate_minhash(&report_stack, ug, 21, 1000, n_thread);
    fprintf(stderr, "[M::%s] deduplicated rescued circles, used %.2fs\n", __func__, Get_T()-time);

    time = Get_T();
    of=0, l=0, ofidx=0;
    if (!is_print_found){
        sprintf(fname, "%s.rescue.fa", asm_opt.output_file_name);
        fp = fopen(fname, "w"); assert(fp);
    }
    while (of<report_stack.n){  // TODO: parallel if slow
        l = report_stack.a[of];
        if (is_print_found){
            fprintf(stderr, "[debug::%s] ===circ%d:\n", __func__, ofidx+1);
            for (int i=0; i<l; i++){
                fprintf(stderr, "[debug::%s]   tig%.6d%c\n", __func__, 
                        (int)(report_stack.a[of+1+i]>>1)+1, "+-"[report_stack.a[of+1+i]&1]);
            }
        }else{
            hamt_ug_write_path_to_fasta(ug, report_stack.a+of+1, l, 1, ofidx, fp);
        }
        of = of+1+l;
        ofidx++;
    }
    fprintf(stderr, "[M::%s] wrote deduplicated rescued circles, used %.2fs\n", __func__, Get_T()-time);

    if (fp) fclose(fp);
    free(fname);
    free(labels);
    free(weights);
    free(color);
    free(used);
    free(scc_counts);

    vu32_t *lt = (vu32_t*)malloc(sizeof(vu32_t));
    kv_init(*lt); 
    kv_resize(uint32_t, *lt, 32);
    for (uint32_t i=0; i<report_stack.n; i++){
        if (ug->u.a[report_stack.a[i]>>1].len>=100000) 
            kv_push(uint32_t, *lt, report_stack.a[i]>>1);
    }
    stacku32_destroy(&report_stack);
    return lt;
}


int hamt_ug_delete_unconnected_single_read_contigs(asg_t *sg, ma_ug_t *ug){
    int ret = 0, verbose = 0;
    asg_t *auxsg = ug->g;
    uint32_t nv;
    asg_arc_t *av;
    ma_utg_t *utg;
    int n1,n2;
    for (uint32_t i=0; i<auxsg->n_seq*2; i++){
        if (i&1) continue;
        if (ug->u.a[i>>1].n==1){
            uint32_t v = ug->u.a[i>>1].start>>1;
            if (verbose) fprintf(stderr, "[d::%s] contig %.6d, start read %.*s. checking..\n", 
                                    __func__, (int)(i>>1)+1,
                                    (int)Get_NAME_LENGTH(R_INF, v), Get_NAME(R_INF, v));
            n1 = n2 = 0;
            nv = asg_arc_n(auxsg, i);
            av = asg_arc_a(auxsg, i);
            for (int j=0; j<nv; j++){
                if (verbose)  fprintf(stderr, "[d::%s]  dir1, target tig %.6d del %d\n", __func__,
                                    (int)(av[j].v>>1)+1, (int)av[j].del );
                if (!av[j].del) n1++;
            }
            if (n1!=0) continue;
            nv = asg_arc_n(auxsg, i^1);
            av = asg_arc_a(auxsg, i^1);
            for (int j=0; j<nv; j++){
                if (verbose) fprintf(stderr, "[d::%s]  dir2, target tig %.6d del %d\n", __func__,
                                    (int)(av[j].v>>1)+1, (int)av[j].del );
                if (!av[j].del) n2++;
            }
            if (n2!=0) continue;

            
            // sg->seq[v].del = 1;
            asg_seq_del(sg, v);
            ret++;
            if (verbose) fprintf(stderr, "[d::%s]    dropping read %.*s\n", __func__,
                        (int)Get_NAME_LENGTH(R_INF, v), Get_NAME(R_INF, v));
            
        } 
    }
    return ret;
}




typedef struct {
    stacku32_t *entries;  // persumed entry vertices, to be tested; do NOT have the direction bit
    stacku32_t *exits;  // (similar)
    stacku32_t *paths;  // length=n_threads, an array of stacks to store paths through tangles
    vu128_t *paths_sten;  // stores entry-exit pairs that are actually used,
                              // to be sorted and ensure stable tangle popping order.
    hamt_ba_t *is_handle; // a bit array, 0 not entry/exit, 1 yes
    ma_ug_t *ug;  // will modify
    asg_t *sg; 
    int n_threads, base_label, gc_tangle_max_tig;
}hamt_markt_t;

typedef struct{
    asg_t *sg;  // string graph
    ma_ug_t *ug;
    int n_threads;
    int bufsize;   // collect how many pairs before forwarding to the next step
    int base_label, gc_tangle_max_tig;
    uint32_t starter, offset;
    uint32_t treated;
    hamt_ba_t *is_handle;
    int debug_step;
}hamt_marknpop_pl_t;

typedef struct{
    hamt_marknpop_pl_t *p;
    hamt_markt_t *dd;
}hamt_marknpop_st_t;


/**
 * @brief Record a walk through a known SESE region, will not modify the graph.
 *        Call the treatwalk method to modify graph.
 *        It is a DFS that greedily pick ndoes with highest coverage at each
 *        step. When seeing a loop, the search will step back, and while 
 *        retreating from a node, we mark all child nodes of this node 
 *        as unvisited. 
 * @note Will not try to early terminate when seeing handle. Let the caller
 *       deal with nesting.
 * @par source forward direction
 * @par sink closing backward, i.e. source->...<-sink.
 * @par pathlength Output, will store lenght of path (counting source and sink)
 * @return 0 if found and recorded path, 1 otherwise.
*/
int hamt_ug_resolveTangles_threaded_callback_markwalk(asg_t *sg, ma_ug_t *ug, 
                                                uint32_t source, uint32_t sink, 
                                                stacku32_t *pathholder,
                                                int *pathlength,
                                                int base_label, int alt_label){
    int verbose = 0;
    uint32_t vu, nv, wu, uu, nu;
    asg_arc_t *av;
    asg_arc_t *au;
    asg_t *auxsg = ug->g;
    int ret = 0;
    *pathlength = 0;
    char debug_name[50];
    sprintf(debug_name, "[%.6d(%c)-%.6d(%c)]", (int)(source>>1)+1, "+-"[source&1],
                        (int)(sink>>1)+1, "+-"[sink&1]);

    // special case, early termination: self circle (shortcut or 1 unitig circle)
    nv = asg_arc_n(auxsg, source);
    av = asg_arc_a(auxsg, source);
    { int cnt=0, self=0; 
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            if (av[i].v==source){
                self = 1;
            }
            cnt++;
        }
        if (self==1 && cnt==1){
            if (verbose) {fprintf(stderr, "[debug::%s] %s special case failure\n", __func__, debug_name);}
            return 1;
        }
    }

    hamt_ba_t *color = hamt_ba_t_init(auxsg->n_seq);

    // a buffer for greedily choosing the path with most coverage
    stacku64_t for_sorting;
    stacku64_init(&for_sorting);
    uint64_t packed_value;
    hamt_ba_t *visited = hamt_ba_t_init(auxsg->n_seq*2);

    stacku32_t s;
    stacku32_init(&s);
    int is_exhausted, is_loop, san=0;
    int is_hinted_circular = ((sink>>1) == (source>>1));
    int is_reached_sink = 0;

    stacku32_push(&s, source);
    hamt_ba_t_write(visited, 1, source);
    hamt_ba_t_write(visited, 1, source^1);
    hamt_ba_t_write(visited, 1, sink);
    hamt_ba_t_write(visited, 1, sink^1);

    while (stacku32_pop(&s, &vu)){
        if (verbose) {fprintf(stderr, "[debug::%s] %s at utg%.6d dir %d\n", 
                                __func__, debug_name, (int)(vu>>1)+1, (int)(vu&1));}
        if (vu==(sink^1)){
            if (is_hinted_circular && (s.n==0)){  // special case: 1st pop of a hinted-circular subgraph, proceed to check child nodes
                ;
            }else{  // regular case: found end, end the walk
                stacku32_push(&s, vu);  // put back
                if (verbose) {fprintf(stderr, "[debug::%s] %s reached sink\n", __func__, debug_name);}
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
        //   (the sorting)
        stacku64_reset(&for_sorting);
        for (int i=0; i<nv; i++){
            if (av[i].del) {continue;}
            if (auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
            packed_value = ((uint64_t)ug->utg_coverage[av[i].v>>1])<<32 | (uint64_t)av[i].v;
            stacku64_push(&for_sorting, packed_value);
            if (verbose){
                fprintf(stderr, "[debug::%s] %s (add to sorting buffer) vertex %.6d , cov %.d\n", 
                    __func__, debug_name,
                    (int)(av[i].v>>1)+1, (int)ug->utg_coverage[av[i].v>>1]);
            }
        }
        radix_sort_ovhamt64(for_sorting.a, for_sorting.a+for_sorting.n);
        stacku64_invert_buffer(&for_sorting);

        // topology check
        for (int i=0; i<for_sorting.n; i++){  // greedily best
            wu = (uint32_t)for_sorting.a[i];
            if (hamt_ba_t_get(visited, wu)) {
                if (wu==(sink^1)){
                    if (verbose) {fprintf(stderr, "[debug::%s] %s reached sink(a)\n", __func__, debug_name);}
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    is_reached_sink = 1;
                    break;
                }
                if (verbose) fprintf(stderr, "[debug::%s] %s   (utg%.6d already visited)\n", __func__, debug_name, (int)(wu>>1)+1);
                continue;
            }
            if (hamt_ba_t_get(visited, wu^1)) {
                if (wu!=(sink^1)){
                    is_loop = 1;
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s] %s reached sink(b)\n", __func__, debug_name);}
                    is_reached_sink = 1;
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    break;
                }
            }
            if (verbose) fprintf(stderr, "[debug::%s] %s    not exhausted: %.6d dir %d\n", __func__, debug_name,
                                (int)(wu>>1)+1, (int)(wu&1));
            is_exhausted = 0;
            break;
        }

        if (is_reached_sink){break;}
        if (is_exhausted && is_loop){  // remove current node
            if (verbose) {fprintf(stderr, "[debug::%s] %s    found looping(a) at utg%.6d, remove %.6d\n", 
                                    __func__, debug_name, (int)(wu>>1)+1, (int)(vu>>1)+1);}
                continue;
        }
        if (!is_exhausted){  // step
            if (hamt_ba_t_get(visited, wu^1)) {  // check if it's a loop
                if (wu!=(sink^1)){
                    if (verbose) {fprintf(stderr, "[debug::%s] %s    found looping(b) at utg%.6d, remove %.6d\n", 
                                    __func__, debug_name, (int)(wu>>1)+1, (int)(vu>>1)+1);}
                    continue;
                }else{
                    if (verbose) {fprintf(stderr, "[debug::%s] %s reached sink(c)\n", __func__, debug_name);}
                    is_reached_sink = 1;
                    stacku32_push(&s, vu);  // put back current round's handle
                    stacku32_push(&s, wu);
                    break;
                }
            }else{
                if (verbose) {fprintf(stderr, "[debug::%s] %s     push %.6d dir %d\n", 
                                        __func__, debug_name, (int)(wu>>1)+1, (int)(wu&1));}
                stacku32_push(&s, vu);  // put back
                stacku32_push(&s, wu);
                hamt_ba_t_write(visited, 1, wu);
            }
        }else{  // remove current node
            // mark all children as unvisited
            for (int i=0; i<nv; i++){
                if (av[i].del) {continue;}
                if (auxsg->seq_vis[av[i].v>>1]!=base_label){continue;}
                hamt_ba_t_write(visited, 0, av[i].v);
                if (verbose) {fprintf(stderr, "[debug::%s] %s       removed child %.6d\n", 
                                               __func__, debug_name, (int)(av[i].v>>1)+1);}
            }
            if (verbose) {fprintf(stderr, "[debug::%s] %s     removed %.6d\n", __func__, debug_name, (int)(vu>>1)+1);}
            ;
        }
        san++;
        if (san>1000){
            fprintf(stderr, "[W::%s] too many steps, abort this try (source was utg%.6d, sink was utg%.6d)\n", __func__, (int)(source>>1)+1, (int)(sink>>1)+1);
            ret = 1;
            goto finish;
        }
    }


    if (s.n==0){
        if (verbose) {fprintf(stderr, "[E::%s] %s can't find a path\n", __func__, debug_name);}
        ret = 1;
        goto finish;
    }else{
        if (verbose){
            if (verbose) {fprintf(stderr, "[debug::%s] %s the path (%d nodes):\n", __func__, debug_name, s.n);}
            for (int i=0; i<s.n; i++){
                fprintf(stderr, " %s       utg%.6d\n", debug_name, (int)(s.a[i]>>1)+1);
            }
        }
        for (int i=0; i<s.n; i++){
            stacku32_push(pathholder, s.a[i]);
        }
        stacku32_push(pathholder, UINT32_MAX);  // marks end of path
        *pathlength = s.n;   
    }

finish:
    stacku32_destroy(&s);
    hamt_ba_t_destroy(visited);
    stacku64_destroy(&for_sorting);
    hamt_ba_t_destroy(color);
    return ret;
}


/**
 * @brief Given a walk through a SESE region, drop vertices and edges beside
 *        those in the walk.
 * @par s (Will be consumed.) A stack storing the path to retain; first and last
 *        is the entry and the exit node. 
*/
void hamt_ug_resolveTangles_threaded_callback_treatwalk(asg_t *sg, ma_ug_t *ug, 
                                                uint32_t source, uint32_t sink, 
                                                stacku32_t *s,
                                                int base_label, int alt_label){
    uint32_t vu, wu, nv, uu, nu;
    asg_arc_t *av, *au;
    asg_t *auxsg = ug->g;

    stacku32_pop(s, &vu); 

    // recover the stored path
    while (stacku32_pop(s, &wu)){
        if (wu!=source && wu!=sink){  // update blacklist here
            ;
        }
        
        hamt_ug_arc_del_selfcircle(sg, ug, wu, base_label);
        
        hamt_ug_arc_del(sg, ug, wu, vu, 0);
        nv = asg_arc_n(auxsg, wu);
        av = asg_arc_a(auxsg, wu);
        for (uint32_t i=0; i<nv; i++){
            if (av[i].v!=vu) {
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

}

/**
 * @brief Callback for kt_for. Given a pair of handles, 1) determine if 
 * they form a single-entry single-exist region, 2) record the path (will not
 * modify).
 * @par data struct
 * @par i # of entry in the input
 * @par tid # of thread
*/
static void hamt_ug_resolveTangles_threaded_callback_step1(void *data, long i, int tid){
    int verbose = 0;
    hamt_markt_t *dd = (hamt_markt_t*) data;
    ma_ug_t *ug = dd->ug;
    asg_t *auxsg = dd->ug->g;

    uint32_t vu;
    uint32_t wu;

    uint64_t pair;
    int max_length, max_coverage, found;
    uint32_t source, sink;



    for (uint32_t dir =0; dir<4; dir++){  // try all direction combinations
        // special case: when entry==exit, 
        //               do not allow self-loop, and check only one direction.
        if (dd->entries->a[i]==dd->exits->a[i]){
            if (dir!=1) continue;  // source->...<sink
        }

        // if (dd->entries->a[i]+1==993 || dd->entries->a[i]==95) verbose = 2;  // 990

        vu = (dd->entries->a[i] <<1) | (dir&1);
        wu = (dd->exits->a[i] <<1) | (dir>>1);

        found = 0;
        if (verbose>1){
            fprintf(stderr, "[debug::%s] vu %.6d (%d) wu %.6d (%d) | i=%d, tid=%d\n", __func__, 
                    (int)(vu>>1)+1, (int)(vu&1), (int)(wu>>1)+1, (int)(wu&1), 
                    (int)i, (int)tid);
        }

        // determine if it's a entry-exit pair
        max_length = 0;
        max_coverage = 0;
        if (hamt_ug_check_localSourceSinkPair(ug, vu, wu, &max_length, &max_coverage, 
                        dd->base_label, dd->gc_tangle_max_tig)){
            source =  vu;
            sink =  wu;
            found = 1;
        }

        // yes, so find a path through the region and save it.
        if (found){  
            int pathlen;
            // note/TODO/bug: this is an arbitrary waterproof, trying to 
            // prevent weird cases (since we searched both directions 
            // during traversal to tolerate inversions/loops)
            int handlemaxlen = ug->u.a[source>>1].len > ug->u.a[sink>>1].len? ug->u.a[source>>1].len : ug->u.a[sink>>1].len;
            if (((ug->u.a[source>>1].len<50000) && (ug->u.a[sink>>1].len<50000)) ||
                ( (max_length>500000 && max_length>handlemaxlen)) ||
                (max_coverage > (ug->utg_coverage[source>>1]*2)) ||
                (max_coverage > (ug->utg_coverage[sink>>1]  *2))){
            // if (0){  // (for debug: allow any length and coverage)
                    if (verbose){
                        fprintf(stderr, "[debug::%s] pair %.6d (%d) and %.6d (%d), but failed length or cov. Max l %d, max cov %d (handles' are %d and %d)\n", 
                                            __func__, 
                                            (int)(source>>1)+1, (int)(source&1), (int)(sink>>1)+1, (int)(sink&1),
                                            max_length, max_coverage, (int)(ug->utg_coverage[source>>1]), (int)ug->utg_coverage[sink>>1]);
                    }
            }else{  // success
                if (verbose){
                    fprintf(stderr, "[debug::%s] found pair %.6d (%d) and %.6d (%d)\n", 
                                            __func__, (int)(source>>1)+1, (int)(source&1), (int)(sink>>1)+1, (int)(sink&1)
                    );
                    }
                // Find a path and store it in thread-specific buffers, 
                // separated by 0xffffffff.
                uint64_t tmpstart = dd->paths[tid].n;  // marks the start of a segment in the path buffer
                hamt_ug_resolveTangles_threaded_callback_markwalk(dd->sg, ug, 
                            source, sink, &dd->paths[tid], &pathlen, 0, 1);
                if (pathlen>0){
                    hamt_ba_t_write(dd->is_handle, 1, source);
                    hamt_ba_t_write(dd->is_handle, 1, sink);
                    kv_push(hamt_u128_t, dd->paths_sten[tid], 
                            ((hamt_u128_t){.x1=(((uint64_t)source)<<32 | sink) , 
                             .x2=(tmpstart <<32 | tid)
                             }
                            )
                    );
                }
                // if (verbose){
                //     fprintf(stderr, "[debug::%s] tid=%d, path buffer length now %d\n", 
                //             __func__, tid, dd->paths[tid].n);
                // }
                
            }
        }

    }
}


/**
 * @brief Callback for kt_pipeline; multithread tangle marking and popping.
 * @note The pipeline is mainly for reducing peak 
 *       memory usage which is quadratic to the number of entry-exit pairs. 
 * @note Popping when only given the entry-exit pair is not fast enough on 
 * single thread, need to parallel it along with the marking (which is 
 * O(|V|**2) for each disconnected subgraph).
 
*/
static void *hamt_ug_resolveTangles_threaded_callback(void *data, 
                                                        int step, void *in){
    hamt_marknpop_pl_t *p = (hamt_marknpop_pl_t*)data;
    double startTime = Get_T();
    int verbose = 0;
    if (step==0){  // mark pairs, also record the path to retain
        int ret;
        // init step data
        hamt_marknpop_st_t *s = (hamt_marknpop_st_t*)calloc(1, sizeof(*s));
        s->p = p;
        s->dd = (hamt_markt_t*) calloc(1, sizeof(hamt_markt_t));
        s->dd->sg = p->sg;
        s->dd->ug = p->ug;
        s->dd->n_threads = p->n_threads;
        s->dd->base_label = p->base_label;
        s->dd->gc_tangle_max_tig = p->gc_tangle_max_tig;
        s->dd->is_handle = p->is_handle;
        s->dd->paths = (stacku32_t*)calloc(p->n_threads, sizeof(stacku32_t));
        s->dd->paths_sten = (vu128_t*)calloc(p->n_threads, sizeof(vu128_t));
        for (int i=0; i<p->n_threads; i++){
            stacku32_init(&s->dd->paths[i]);
            kv_init(s->dd->paths_sten[i]);
            kv_resize(hamt_u128_t, s->dd->paths_sten[i], 16);
        }
        s->dd->entries = (stacku32_t*)calloc(1, sizeof(stacku32_t));
        s->dd->exits = (stacku32_t*)calloc(1, sizeof(stacku32_t));
        stacku32_init(s->dd->entries);
        stacku32_init(s->dd->exits);

        // fill in handles to be tested
        {
            uint32_t entry;
            uint32_t exit = p->starter + p->offset;
            for (entry=p->starter; entry<p->ug->g->n_seq*2; entry+=2, p->starter+=2){
                while (exit<p->ug->g->n_seq*2){
                    if (p->ug->u.a[entry>>1].subg_label==p->ug->u.a[exit>>1].subg_label){  // &&
                        stacku32_push(s->dd->entries, entry>>1);
                        stacku32_push(s->dd->exits, exit>>1);
                        // if (p->debug_step==1) fprintf(stderr, "GO push %.6d %.6d\n", (int)(entry>>1)+1, (int)(exit>>1)+1);
                    }else{
                        ;
                        // if (p->debug_step==0) fprintf(stderr, "WUT skip push %.6d %.6d | entry label %d, exit label %d, isdtip %d %d\n", 
                        //     (int)(entry>>1)+1, (int)(exit>>1)+1,
                        //     (int)p->ug->u.a[entry>>1].subg_label, (int)p->ug->u.a[exit>>1].subg_label,
                        //     hamt_asgarc_util_isDanglingTip(p->ug->g, entry, 0, 0, p->base_label),
                        //     hamt_asgarc_util_isDanglingTip(p->ug->g, exit, 0, 0, p->base_label)
                        //     );
                    }
                    exit+=2;
                    p->offset+=2;
                    if (s->dd->entries->n==p->bufsize) break;
                }
                if (s->dd->entries->n==p->bufsize) break;
                exit = entry+2;
                p->offset = 0;
            }
            if (verbose) fprintf(stderr, "[debug::%s] step0, dd length %d, new entry starts at %.6d , exit at %.6d\n", 
                __func__, s->dd->entries->n,  (int)(entry>>1)+1, (int)(exit>>1)+1);
        }
        return s;
    }else if (step==1){
        hamt_marknpop_st_t *s = (hamt_marknpop_st_t*) in;
        hamt_marknpop_pl_t *p = s->p;
        if (s->dd->entries->n>0){  // find SESE and figure out path
            kt_for(p->n_threads-2<1? 1:p->n_threads-2, hamt_ug_resolveTangles_threaded_callback_step1,
                        s->dd, s->dd->entries->n);
            return s;
        }else{  // end of stream
            stacku32_destroy(s->dd->entries);
            stacku32_destroy(s->dd->exits);
            free(s->dd->exits);
            free(s->dd->entries);
            for (int i=0 ;i<p->n_threads; i++){
                stacku32_destroy(&s->dd->paths[i]);
                kv_destroy(s->dd->paths_sten[i]);
            }
            free(s->dd->paths);
            free(s->dd->paths_sten);
            free(s->dd);
            free(s);
        }

    }else if (step==2){  // flip 
        hamt_marknpop_st_t *s = (hamt_marknpop_st_t*) in;
        hamt_marknpop_pl_t *p = s->p;
        // To ensure stable tangle popping order, sort entry-exit pairs here.
        // Order is arbitrary (by contig IDs).
        vu128_t orderbuf;
        kv_init(orderbuf);
        kv_resize(hamt_u128_t, orderbuf, 16);
        for (int tid=0; tid<p->n_threads; tid++){
            for (int i=0; i<s->dd->paths_sten[tid].n; i++){
                kv_push(hamt_u128_t, orderbuf, s->dd->paths_sten[tid].a[i]);
            }
        }
        qsort(orderbuf.a, orderbuf.n, sizeof(hamt_u128_t), compare_u128_t_use64);
        
        // SESE regions can only contain or exclude
        //  each other, but they never partially overlap. 
        // Here, do not treat SESE that contains any other SESE.
        stacku32_t tmp;
        stacku32_init(&tmp);
        stacku32_t *stack;
        stacku32_t sancheck_contain;
        stacku32_init(&sancheck_contain);
        if (verbose) {
            fprintf(stderr, "[dbg::%s] oderbuf length=%d\n", __func__, (int)orderbuf.n);
            for (int i_buf=0; i_buf<orderbuf.n; i_buf++){
                fprintf(stderr, "[dbg::%s]   %d : %d / %d\n", __func__, 
                        i_buf,
                        (int)((uint32_t)(orderbuf.a[i_buf].x1>>32)), 
                        (int)((uint32_t)orderbuf.a[i_buf].x1)
                        );
            }
        }
        for (int i_buf=0; i_buf<orderbuf.n; i_buf++){
            uint32_t tid =  (uint32_t)orderbuf.a[i_buf].x2;
            uint32_t buf_st = orderbuf.a[i_buf].x2>>32;

            int has_contain = 0;
            stack = &s->dd->paths[tid];
            for (int i=buf_st; i<stack->n; i++){
                if (stack->a[i]==UINT32_MAX){
                    if (hamt_ba_t_get(s->dd->is_handle, stack->a[i-1])){
                        has_contain-=1;
                    }
                    if (has_contain==0 && tmp.n>0){
                        if (verbose){
                            fprintf(stderr, "[debug::%s] treated source %.6d (%d) sink %.6d (%d), size=%d, $has_contain=%d\n", 
                                    __func__, (int)(tmp.a[0]>>1)+1, (int)(tmp.a[0]&1),
                                    (int)(tmp.a[tmp.n-1]>>1)+1, (int)(tmp.a[tmp.n-1]&1),
                                    tmp.n,
                                    has_contain
                                    );
                            for (int tmpi=0; tmpi<tmp.n; tmpi++){
                                fprintf(stderr, "[debug::%s]\t+\t%.6d\n", __func__, 
                                            (int)(tmp.a[tmpi]>>1)+1);
                            }
                        }
                        // treat
                        hamt_ug_resolveTangles_threaded_callback_treatwalk(
                            p->sg, p->ug, 
                            tmp.a[0], tmp.a[tmp.n-1], &tmp, 
                            p->base_label, p->base_label^1);
                        p->treated++;
                        
                    }else{
                        if (verbose){
                            fprintf(stderr, "[debug::%s] ignore source %.6d (%d) sink %.6d (%d); size=%d, $has_contain=%d\n", 
                                    __func__, (int)(tmp.a[0]>>1)+1, (int)(tmp.a[0]&1),
                                    (int)(tmp.a[tmp.n-1]>>1)+1, (int)(tmp.a[tmp.n-1]&1),
                                    tmp.n, 
                                    has_contain
                                    );
                            for (int tmpi=0; tmpi<sancheck_contain.n; tmpi++){
                                fprintf(stderr, "[debug::%s]\tm\t%.6d\n", __func__, 
                                            (int)(sancheck_contain.a[tmpi]>>1)+1);
                            }
                        }
                    }
                    has_contain = 0;
                    stacku32_reset(&tmp);
                    stacku32_reset(&sancheck_contain);
                    break;
                }else{
                    if (hamt_ba_t_get(s->dd->is_handle, stack->a[i]) && 
                        tmp.n!=0){
                        has_contain+=1;
                        stacku32_push(&sancheck_contain, stack->a[i]);
                    }
                    stacku32_push(&tmp, stack->a[i]);
                }
            }
        }

        kv_destroy(orderbuf);
        stacku32_destroy(&sancheck_contain);
        stacku32_destroy(s->dd->entries);
        free(s->dd->entries);
        stacku32_destroy(s->dd->exits);
        free(s->dd->exits);
        for (int i=0 ;i<p->n_threads; i++){
            stacku32_destroy(&s->dd->paths[i]);
            kv_destroy(s->dd->paths_sten[i]);
        }
        kv_destroy(tmp);
        free(s->dd->paths);
        free(s->dd->paths_sten);
        free(s->dd);
        free(s);

    }else{
        fprintf(stderr, "[E::%s] does not have such step, check caller code.\n", __func__);
        exit(1);
    }
    return 0;
}



/**
 * @brief Multithreaded tangle popping.
 * @note Currently O(n**2) for marking; reason is that SESE provides O(|V|)
 * but it requires ordinary SCC, which on the aux graph does not tolerate 
 * inversions on the bidirected graph. In practice, if there is a huge tangle
 * in the assembly graph, SCC regions are usually quite limited. 
 * @par sg string graph
 * @par ug unitig/contig graph
 * @par n_threads # of threads
 * @par base_label 0 to work with primary graph
 * @return # spots treated; will print stats to stderr.
*/
int hamt_ug_resolveTangles_threaded(asg_t *sg, ma_ug_t *ug, int n_threads, int base_label, int debug_step){
    asg_t *auxsg = ug->g;
    double startTime = Get_T();

    // prepare input
    int max_label = hamt_ug_util_BFS_markSubgraph(ug, base_label);
    hamt_marknpop_pl_t *p = (hamt_marknpop_pl_t*)calloc(1, sizeof(hamt_marknpop_pl_t));
    p->sg = sg;
    p->ug = ug;
    p->n_threads = n_threads;
    p->bufsize = 48*4096;  // must not be variable
    p->base_label = base_label;
    p->gc_tangle_max_tig = asm_opt.gc_tangle_max_tig;
    p->is_handle = hamt_ba_t_init(auxsg->n_seq*2);
    p->debug_step = debug_step;

    // mark and pop
    kt_pipeline(3, hamt_ug_resolveTangles_threaded_callback, 
                p, 3);

    // update graph
    if (p->treated>0){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
        hamt_ug_cleanup_arc_by_labels(sg, ug);
    }


    // final
    fprintf(stderr, "[M::%s]   cleaning used %.1fs\n", __func__, Get_T()-startTime);
    fprintf(stderr, "[M::%s] treated %d spots\n", __func__, p->treated);
    int ret = p->treated;

    // clean up
    hamt_ba_t_destroy(p->is_handle);
    free(p);
    return ret;
}


/**
 * @brief Disconnect two long contigs if their coverages are very different.
 * @return # of treatments.
*/
int hamt_ug_disconnect_long_contig_pairs_by_cov(asg_t *sg, ma_ug_t *ug){
    int lt = 200000;   // length threshold
    int lt_tip = 100000;  // less stringent for tips
    int ct_maxmin = 30;  // coverage threshold, minimum of the longer of the two
    int ct_maxmin_tip = 10;   // less stringent if one is a tip
    int ct_minratio = 3;  // coverage ratio threshold,  not inclusive

    int n_treated = 0;
    double T = Get_T();
    uint32_t vu, wu, lv, lw, cv, cw;  // id, id, length, length, coverage, coverage
    int vu_istip, wu_istip;
    uint32_t nv;
    asg_arc_t *av;
    asg_t *auxsg = ug->g;

    int tlt, tct;  // tmp lt, tmp ct-maxmin
    for (vu=0; vu<auxsg->n_seq*2; vu++){
        lv = ug->u.a[vu>>1].len;
        vu_istip = hamt_asgarc_util_isTip(auxsg, vu, 0,0, -1);
        
        nv = asg_arc_n(auxsg, vu);
        av = asg_arc_a(auxsg, vu);
        for (int i=0; i<nv; i++){
            if (av[i].del) continue;
            wu = av[i].v;
            lw = ug->u.a[wu>>1].len;
            wu_istip = hamt_asgarc_util_isTip(auxsg, wu, 0,0, -1);

            cv = ug->utg_coverage[vu>>1];
            cw = ug->utg_coverage[wu>>1];

            // check length and coverage
            if (vu_istip){
                tlt = lt_tip;  tct = ct_maxmin_tip;
            }else{
                tlt = lt;  tct = ct_maxmin;
            }
            if (lv<tlt || cv<tct) continue;
            if (wu_istip){
                tlt = lt_tip;  tct = ct_maxmin_tip;
            }else{
                tlt = lt;  tct = ct_maxmin;
            }
            if (lw<tlt || cw<tct) continue;
            //   (coverage ratio)
            float ratio;
            if (cv>cw) ratio = ((float)cv)/cw;
            else ratio = ((float)cw)/cv;
            if (ratio<=ct_minratio) continue;

            // delete the arc
            hamt_ug_arc_del_between(sg, ug, vu, wu, 1, 0); 
            n_treated++;
        }
    }

    // update graph
    if (n_treated>0){
        asg_cleanup(sg);
        asg_cleanup(auxsg);
    }

    fprintf(stderr, "[M::%s] treated %d, used %.1fs.\n", __func__, n_treated, 
            Get_T()-T);
    return n_treated;
}

/**
 * @par buf At least is 512 long.
 * @func Collect canonical 5-mer profile of the sequence.
 * 
*/
void hamt_5NF_profile_from_seq_core(char *seq, uint32_t seq_l, double *buf){
    memset(buf, 0, sizeof(double)*512);
    uint16_t mer=0;
    uint16_t mer_rev = 0;
    uint32_t mer_l=0;
    uint16_t c;
    uint16_t mer_mask = (uint16_t)(1<<10)-1;

    uint32_t tot = 0;
    for (uint32_t i=0; i<seq_l-5+1; i++){
        c = (uint16_t)seq_nt4_table[(int)seq[i]];
        if (c<4){
            mer <<= 2;
            mer |= c;
            mer &= mer_mask;
            //mer_rev>>=2;
            //mer_rev = c<<8; 
            //mer_rev &=mer_mask;
            //fprintf(stderr, "dbg %d %d\n", (int)mer, (int)mer_rev);
            if (mer_l==4){
                buf[index5NF[mer]]++;
                tot++;
            }else mer_l++;
        }else{
            mer = mer_l = 0;
        }
    }
    for (int i=0; i<512; i++){
        buf[i] /= tot;
    }
}
typedef struct{
    ma_ug_t *ug;
    uint64_t *tigIDs;  // will only use the right 32bits
                       // (is this wrong if different endianess??)
    uint32_t tigIDs_l;
    double *mat;
}hamt_get5nf_t;
static void hamt_5NFwithcov_profile_from_seq_callback(void *data, long jobID, int threadID){
    hamt_get5nf_t *s = (hamt_get5nf_t*)data;
    uint32_t v = (uint32_t)s->tigIDs[jobID];
    char *seq = s->ug->u.a[v].s;
    uint32_t seq_l = s->ug->u.a[v].len;
    double *buf = s->mat+(513*jobID);
    buf[0] = log10((double)s->ug->utg_coverage[v]);  // note: this needs to be collect. hamt ug regen always does it.
    hamt_5NF_profile_from_seq_core(seq, seq_l, buf+1);
}
void hamt_5NF_profile_gen(ma_ug_t *ug, uint64_t *candidates, uint32_t candidates_l, 
        int n_threads, double **mat){
    hamt_get5nf_t s;
    s.ug = ug;
    s.tigIDs = candidates;
    s.tigIDs_l = candidates_l;
    s.mat = *mat;
    kt_for(n_threads, hamt_5NFwithcov_profile_from_seq_callback, &s, candidates_l);

    // normalize
    znorm2D(mat, (int)candidates_l, 513);

}

int compare_ddp_t(const void *a_, const void *b_){
    ddp_t *a = (ddp_t*)a_;
    ddp_t *b = (ddp_t*)b_;
    if (a->is_optimal && !b->is_optimal) return -1;
    else if (!a->is_optimal && b->is_optimal) return 1;

    if (a->d1 > b->d1) return 1;
    else if (a->d1 < b->d1) return -1;
    else{
        if (a->d1 < b->d1) return 1;  // reverse sort
        else if (a->d1 > b->d1) return -1;
    }
    return 0;
}
void hamt_simple_binning_pick_neighbors(FILE *fp, char *bin_dir_prefix, ma_ug_t *ug, 
                                        vu64_t *candidates, double *emb,
                                        double radius1, double radius2, 
                                        int *counter_tot, int*counter_mul){  

    vddp_t itvls;
    kv_init(itvls);
    kv_resize(ddp_t, itvls, 16);

    // (upper 32 bits of candidates were used for sorting; 
    // contig ID and are always ui32 in this function,
    // so we just rely on unsigned truncation.)

    // blocking candidates:
    //    - block: wrt idx of the candidates buffer, is used for all iterations 
    //    - the .is_ignored property will be reset for each clustering step,
    //      as we re-collect angle intervals for each new seed.
    //    (will need to check && update both)
    int verbose = 0;

    uint8_t *block = (uint8_t*)calloc(candidates->n, 1);  // can recycle upper 32bits of candidates instead, but doesn't matter
    assert(block);

    int I = 0; // binID

    uint32_t seedID, otherID;
    double seedx, seedy, x, y, r, t;
    double boundary=radius2*2;
    for (int32_t i=candidates->n-1; i>=0; i--){  // seed from the longest candidate contigs
        if (verbose) fprintf(stderr, "[dbg::%s] seed ctg%.6d, i is %d\n", 
                __func__, (int)((uint32_t)candidates->a[i])+1, (int)i);
        if (block[i]) continue;
        itvls.n = 0;

        seedID = candidates->a[i];
        seedx = emb[2*i];
        seedy = emb[2*i+1];

        for (int32_t j=i-1; i>=1 && j>=0; j--){
            if (block[j]) continue;

            // to polar
            x = emb[2*j] - seedx;  // shift origin point wrt seed
            y = emb[2*j+1] - seedy;
            r = sqrt(pow(x, 2) + pow(y, 2));

            if (r<boundary){
                t = atan(y/(x+0.000000001));
                if (x<0) t+=HAMT_PI;
                else if (x>=0 && y<0) t+=HAMT_TWOPI;

                // get range of available angle
                double theta = acos(1-pow(r, 2)/(2*pow(radius2, 2)));
                double relative_angle = (HAMT_PI-theta)/2;
                double st = t-relative_angle;
                double en = t+relative_angle;
                kv_push(ddp_t, itvls, ((ddp_t){
                                        .d1=st, 
                                        .d2=en, .i=(uint32_t)j,
                                        .is_optimal=r<radius1?(uint8_t)1:(uint8_t)0, 
                                        .is_ignored=0}));
                if (st<0){  // trying to patch where it wraps around. maybe have bug
                    kv_push(ddp_t, itvls, ((ddp_t){
                                        .d1=st+HAMT_TWOPI, 
                                        .d2=en+HAMT_TWOPI, .i=(uint32_t)j,
                                        .is_optimal=r<radius1?(uint8_t)1:(uint8_t)0, 
                                        .is_ignored=0}));
                }
            }
        }
        if (verbose) fprintf(stderr, "[dbg::%s]    itvls len=%d\n",__func__, (int)itvls.n);

        // sort, try optimal solution, then remove containment and duplicates and see if suboptimal exists
        qsort(itvls.a, itvls.n, sizeof(ddp_t), compare_ddp_t);
        // (check if there's optimal choice, if so, we are done)
        {
            int n = 0;
            for (int tmpi=0; tmpi<itvls.n; tmpi++){
                if (block[itvls.a[tmpi].i]) continue;
                if (itvls.a[tmpi].is_optimal) n++;
                else break;
            }
            if (n>0){  // write content and update blacklist, then goes to the next seed
                if (verbose) fprintf(stderr, "[dbg::%s]     optimal\n", __func__);
                if (counter_tot) (*counter_tot)++;
                fprintf(fp, "bin%d.opt\ts%d.ctg%.6d%c", I, (int)ug->u.a[seedID].subg_label, 
                            seedID+1, "lc"[ug->u.a[seedID].circ]);
                block[i] = 1;
                FILE *fp_bin = 0;
                if (bin_dir_prefix){
                    char *bin_fn = (char*)malloc(strlen(bin_dir_prefix)+50);
                    sprintf(bin_fn, "%s/bin%d.opt.fa", bin_dir_prefix, I);
                    fp_bin = fopen(bin_fn, "w");
                    assert(fp_bin);
                    free(bin_fn);
                    uint32_t v = seedID;
                    fprintf(fp_bin, ">s%d.ctg%.6d%c\n", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                    fprintf(fp_bin, "%.*s\n", (int)ug->u.a[v].len, ug->u.a[v].s);
                }

                uint32_t i_can;
                for (i_can=0; i_can<n; i_can++){
                    if (block[itvls.a[i_can].i]) continue;
                    uint32_t v = candidates->a[itvls.a[i_can].i];
                    fprintf(fp, "\ts%d.ctg%.6d%c", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                    if (fp_bin){
                        fprintf(fp_bin, ">s%d.ctg%.6d%c\n", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                        fprintf(fp_bin, "%.*s\n", (int)ug->u.a[v].len, ug->u.a[v].s);
                    }
                    block[itvls.a[i_can].i] = 1;
                }
                fprintf(fp, "\n");
                if (fp_bin) fclose(fp_bin);
                I++;
                if (i_can>1)
                    if (counter_mul) (*counter_mul)++;  // neighbor(s) + seed itself
                continue;
            }
        }
        // (If we reach here, no candidates were available in seed's vicinity.
        //  We will try finding a circle that intersects seed and contain all 
        //  nearby candidates; if there are more than one circle or candidates
        //  are within the neighborhood cannot be covered by only one such circle, 
        //  we deny all solutions.)
        // (remove interval containments in the "available angles" intervals)
        for (int iseg=0; iseg<itvls.n; iseg++){
            if (block[itvls.a[iseg].i]) continue;
            if (itvls.a[iseg].is_ignored) continue;
            int jseg = iseg+1;
            double end = itvls.a[iseg].d2;
            itvls.a[iseg].weight = 1;
            while (jseg<itvls.n && itvls.a[jseg].d1<end){
                if (block[itvls.a[jseg].i] || itvls.a[jseg].is_ignored) {jseg++; continue;}
                if (itvls.a[jseg].d2<end){
                    itvls.a[jseg].is_ignored = 1;
                    itvls.a[iseg].weight++;
                }
                jseg++;
            }
        }

        // (slide through sorted intervals to find the highest coverage)
        // (If we allow multiple circles to be selected, we can repurpose 
        //  the .is_optimal property to indicate who is currently
        //  within the window's span. When we see a new high score, push 
        //  to a stack. )
        // The requirement is the suboptimal choice should contain
        //  all neighboring points. Therefore here we just count and see if 
        //  we ever reach the highest score.
        //  Note that this REQUIRES the interval lists was collected from 
        //  valid neighbors (i.e. no bounding box, calculate distance wrt seed). 
        //  Although in practic I guess this doesn't matter, binning mistakes are 
        //  from elsewhere...
        int best = -1;  // known best coverage
        int cov = 0;   // current coverage
        double left;  // left of the current window
        int i_left_window = 0;  // left most window that hasn't moved out yet
        for (int iseg=0; iseg<itvls.n; iseg++){
            if (itvls.a[iseg].is_ignored || block[itvls.a[iseg].i]) continue;
            left = itvls.a[iseg].d1;
            cov += itvls.a[iseg].weight;
            // check if some window on the left has moved out
            while (1){
                if (itvls.a[i_left_window].d2<left){
                    cov -= itvls.a[i_left_window].weight;
                    i_left_window++;
                }else break;
            }
            // check score
            best = best>cov? best : cov;
        }
        if (best>=itvls.n){
            if (verbose) fprintf(stderr, "[dbg::%s]     suboptimal\n", __func__);
            if (counter_tot) (*counter_tot)++;
            fprintf(fp, "bin%d.sub\ts%d.ctg%.6d%c", I, (int)ug->u.a[seedID].subg_label, 
                        seedID+1, "lc"[ug->u.a[seedID].circ]);
            block[i] = 1;

            FILE *fp_bin = 0;
            if (bin_dir_prefix){
                char *bin_fn = (char*)malloc(strlen(bin_dir_prefix)+50);
                sprintf(bin_fn, "%s/bin%d.sub.fa", bin_dir_prefix, I);
                fp_bin = fopen(bin_fn, "w");
                assert(fp_bin);
                free(bin_fn);
                uint32_t v = seedID;
                fprintf(fp_bin, ">s%d.ctg%.6d%c\n", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                fprintf(fp_bin, "%.*s\n", (int)ug->u.a[v].len, ug->u.a[v].s);
            }

            uint32_t i_can=0;
            for (uint32_t i_can=0; i_can<itvls.n; i_can++){
                if (block[itvls.a[i_can].i]) continue;
                uint32_t v = candidates->a[itvls.a[i_can].i];
                fprintf(fp, "\ts%d.ctg%.6d%c", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                if (fp_bin){
                    fprintf(fp_bin, ">s%d.ctg%.6d%c\n", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                    fprintf(fp_bin, "%.*s\n", (int)ug->u.a[v].len, ug->u.a[v].s);
                }
                block[itvls.a[i_can].i] = 1;
            }
            fprintf(fp, "\n");
            if (fp_bin) fclose(fp_bin);
            I++;
            if (i_can>1)
                if (counter_mul) (*counter_mul)++;  // neighbor(s) + seed itself
            continue;
        }

        // If we reach here, put the seed contig alone in a bin if it is long.
        // (if the contig is short, do nothing; bins initiated from other seeds might
        // want to grab it.) 
        if (verbose) fprintf(stderr, "[dbg::%s]     fall-through\n", __func__);
        if ((candidates->a[i]>>32)>500000){
            if (verbose) fprintf(stderr, "[dbg::%s]     write because long (%d)\n", __func__, 
                    (int)(candidates->a[i]>>32));
            if (counter_tot) (*counter_tot)++;
            fprintf(fp, "bin%d.sub\ts%d.ctg%.6d%c\n", I, (int)ug->u.a[seedID].subg_label, 
                        seedID+1, "lc"[ug->u.a[seedID].circ]);

            if (bin_dir_prefix){
                uint32_t v = seedID;
                char *bin_fn = (char*)malloc(strlen(bin_dir_prefix)+50);
                sprintf(bin_fn, "%s/bin%d.thro.fa", bin_dir_prefix, I);
                FILE *fp_bin = fopen(bin_fn, "w");
                assert(fp_bin);
                free(bin_fn);
                fprintf(fp_bin, ">s%d.ctg%.6d%c\n", (int)ug->u.a[v].subg_label, v+1, "lc"[ug->u.a[v].circ]);
                fprintf(fp_bin, "%.*s\n", (int)ug->u.a[v].len, ug->u.a[v].s);
                fclose(fp_bin);
            }
            block[i] = 1;
            I++;
        }
    }    
    free(block);
    kv_destroy(itvls);

}

void hamt_helper_write_tab_joined_darray(FILE *fp, double *a, uint32_t n){
    if (n==0) return;
    fprintf(fp, "%f", a[0]);
    for (int i=1; i<n; i++){
        fprintf(fp, "\t%f", a[i]);
    }
    fprintf(fp, "\n");
}
/**
 * @func Post-assembly, post-circle-resuce binning based on bhtsne.
*/
void hamt_simple_binning(ma_ug_t *ug, vu32_t *blacklist, int n_threads, 
                        char *output_prefix, int write_binning_fasta){
    // Take contigs >=100kb and is not [>=1Mb and circular] and is not used by 
    //  circle rescue step (after deduplication).
    // Make a 2D feature matrix {n_contigs, n_features}, where features are:
    //    - coverage (hifiasm's estimation, tot bases / contig length)
    //    - canonical 5-mer profile (no adjustments to frequencies, 
    //      unlike the n=103 TNF used in various binners, since tSNE seems
    //      indifferent to the difference. Prominent problem lies elsewhere.
    //      Also does not matter too much whether we use 4-mer or 5-mer.) 
    double T = Get_T();
    int verbose = 0;
    
    asg_t *auxsg = ug->g;
    vu64_t nl;  // "neighbor list"
    kv_init(nl); 
    kv_resize(uint64_t, nl, 32);

    // collect candidates
    radix_sort_ovhamt32(blacklist->a, blacklist->a+blacklist->n);
    uint32_t b_p=0;  // blacklist pointer
    for (uint32_t v=0; v<auxsg->n_seq; v++){
        if (ug->u.a[v].len>=100000){
            while (v>blacklist->a[b_p]) 
                b_p++;
            int ban = 0;
            for (uint32_t i=b_p; i<blacklist->n; i++){
                if (blacklist->a[i]==v) {
                    ban = 1;
                    break;
                }
                if (blacklist->a[i]>v) break;
            }
            if (!(ug->u.a[v].len>=1000000 && ug->u.a[v].circ) && !ban){
                kv_push(uint64_t, nl, (((uint64_t)ug->u.a[v].len)<<32)|v );  // for sorting
            }
        }
    }

    // check if skipping tsne
    if (nl.n - 1 < 3 * asm_opt.tsne_perplexity){
        fprintf(stderr, "[M::%s] Too few contigs, skip tsne binning.\n", __func__);
        // cleanup
        kv_destroy(nl);
        // (dummy output files)
        char *fn = (char*)malloc(strlen(output_prefix)+20);
        assert(fn);
        sprintf(fn, "%s.bins.tsv", output_prefix);
        FILE *fp = fopen(fn, "w"); 
        assert(fp);
        fclose(fp);
        if (write_binning_fasta){
            sprintf(fn, "%s.bins", output_prefix);
            mkdir(fn, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);  // 755
        }
        free(fn);
        return;
    }

    fprintf(stderr, "[M::%s] Will try to bin on %d contigs (skipped %d because blacklist).\n", 
            __func__, (int)nl.n, (int)blacklist->n);
    radix_sort_ovhamt64(nl.a, nl.a+nl.n);

    // coverage and 5NF profile, and z-score normalize it
    double *pf = (double*)calloc(513*nl.n, sizeof(double));
    assert(pf);
    hamt_5NF_profile_gen(ug, nl.a, nl.n, n_threads, &pf);

    // call tSNE
    double *emb = ts_fit(nl.n, 513, pf, 2, 0.5, asm_opt.tsne_perplexity, asm_opt.tsne_randomseed);

    // (debug: write embedding)
    if (verbose)
    {
        char *fn = (char*)malloc(strlen(output_prefix)+50);
        assert(fn);

        sprintf(fn, "%s.bins.featuremat", output_prefix);
        FILE *fp = fopen(fn, "w");
        for (int tmpi=0; tmpi<nl.n; tmpi++){
            uint32_t tmpv = (uint32_t)nl.a[tmpi];
            fprintf(fp, "s%d.ctg%.6d%c\t", (int)ug->u.a[tmpv].subg_label, 
                    (int)(tmpv+1), "lc"[ug->u.a[tmpv].circ]
                    );
            hamt_helper_write_tab_joined_darray(fp, pf+tmpi*513, 513);
        }

        fclose(fp);

        sprintf(fn, "%s.bins.embedding", output_prefix);
        fp = fopen(fn, "w");
        for (int tmpi=0; tmpi<nl.n; tmpi++){
            uint32_t tmpv = (uint32_t)nl.a[tmpi];
            fprintf(fp, "s%d.ctg%.6d%c\t", (int)ug->u.a[tmpv].subg_label, 
                    (int)(tmpv+1), "lc"[ug->u.a[tmpv].circ]
                    );
            hamt_helper_write_tab_joined_darray(fp, emb+tmpi*2, 2);
        }
        fclose(fp);
        free(fn);
    }

    // parse the embedding
    free(pf);
    char *fn = (char*)malloc(strlen(output_prefix)+20);
    assert(fn);
    sprintf(fn, "%s.bins.tsv", output_prefix);
    FILE *fp = fopen(fn, "w"); 
    assert(fp);
    if (write_binning_fasta){
        sprintf(fn, "%s.bins", output_prefix);
        mkdir(fn, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);  // 755
    }
    int tot_bins=0, tot_bins_multi=0;
    hamt_simple_binning_pick_neighbors(fp, write_binning_fasta? fn : 0,
                                       ug, &nl, emb, 
                                       asm_opt.tsne_neigh_dist, asm_opt.tsne_neigh_dist*0.8, 
                                       &tot_bins, &tot_bins_multi);
    fclose(fp);
    free(fn);

    // cleanup
    kv_destroy(nl);
    free(emb);
    fprintf(stderr, "[M::%s] Binning used %.2fs. %d bins (%d have more than 1 contig).\n", 
            __func__, Get_T()-T, tot_bins, tot_bins_multi);

}



typedef struct {
    //hamt_ba_t *todel;  // note to self: must not use bit array, will data race;
                         // struct bit field will ub data race 
    uint8_t *todel2;// array to log which read should be deleted
    ma_hit_t_alloc* sources;
    ma_sub_t *coverage_cut; 
    R_to_U *ruIndex;
    int max_hang;
    int min_ovlp;
}hamt_mahitcontadv_t;

static void hit_contained_advance_callback(void *data, 
                                    long jobID,   // is readID
                                    int threadID){  // for kt_for
    hamt_mahitcontadv_t *s = (hamt_mahitcontadv_t*)data;  
    if (s->coverage_cut[jobID].del) return;

    //int dbg_print = 0;
    //if (jobID==67310 || jobID==96791 || jobID==161897) dbg_print = 1;

    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;
    int32_t r;
    uint32_t tn;
    asg_arc_t t;
    //int printed = 0;
    for (long long j=0; j<(long long)s->sources[jobID].length; j++) {
        h = &(s->sources[jobID].buffer[j]);
        //check the corresponding two reads 
        tn = Get_tn(*h);
		sq = &(s->coverage_cut[jobID]);
        st = &(s->coverage_cut[tn]);
        if(sq->del || st->del) continue;  // "may have trio bugs"
        if(h->del) continue;  // "may have trio bugs"
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, 
                    s->max_hang, asm_opt.max_hang_rate, 
                    s->min_ovlp, &t);
        if (r == MA_HT_QCONT){
            //
            // need to consider a rare case when updating R_to_U: 
            // identical reads will be considered as QCONT in both directions..
            // (this was not a problem in the original function due to every 
            // containment encounter triggers .del marking in both direction; 
            // here we collect in one direction and only mark when all are 
            // collected.)
            //
            // A bad solution for now: if tn is smaller than qn and hit is 
            // between identical reads, choose to do nothing.
            int is_special_case = 0;
            if (tn<jobID){
                if (sq->s==0 && st->s==0 && sq->e==st->e &&
                    sq->e==Get_READ_LENGTH(R_INF, jobID) &&
                    st->e==Get_READ_LENGTH(R_INF, tn)){

                    is_special_case = 1;
                    //fprintf(stderr, "[W::%s] saw a pair of identical reads: "
                    //        "%.*s (len %d) and %.*s (len %d) | sancheck: "
                    //        "sq-s %d sq-e %d st-s %d st-e %d \n", __func__, 
                    //        (int)Get_NAME_LENGTH(R_INF, jobID), Get_NAME(R_INF, jobID),
                    //        (int)Get_READ_LENGTH(R_INF, jobID),
                    //        (int)Get_NAME_LENGTH(R_INF, tn), Get_NAME(R_INF, tn),
                    //        (int)Get_READ_LENGTH(R_INF, tn),
                    //        (int)sq->s, (int)sq->e, (int)st->s, (int)st->e
                    //        );

                }
            }
            if (!is_special_case){
                s->todel2[jobID] = 1;//hamt_ba_t_write(s->todel, 1, Get_qn(*h));
                h->del = 1;
                uint32_t tmpleft=jobID, tmpright, tmpisunitig;
                set_R_to_U(s->ruIndex, Get_qn(*h), Get_tn(*h), 0);  // safe in this direction
            }
 
            //if (!printed){
            //    fprintf(stderr, "[dbg::%s] QCONT qn %d tn %d\n", __func__, 
            //            (int)Get_qn(*h), (int)Get_tn(*h));
            //    printed=1;
            //}
        }else if (r == MA_HT_TCONT){
            h->del = 1;  // DO NOTHING about the target: let it happen as QCONT when tn is treated as self
        }
    }
}

int hamt_ma_hit_contained_advance(ma_hit_t_alloc* sources, long long n_read, 
                                  ma_sub_t *coverage_cut, 
                                  R_to_U* ruIndex, int max_hang, int min_ovlp)
{
    // ma_hit_contained_advance improved for speed.
    double startTime0 = Get_T();
    double startTime = Get_T();
    int verbose = 0;

    uint8_t *todel2 = (uint8_t*)calloc(n_read, 1);

    int ret0 = 0;
    int ret = 0;
	int32_t r;
	long long i, j, m;
	asg_arc_t t;
    ma_hit_t *h = NULL;
    ma_sub_t *sq = NULL;
    ma_sub_t *st = NULL;

    hamt_mahitcontadv_t dworker;
    dworker.todel2 = todel2;
    dworker.sources = sources;
    dworker.coverage_cut = coverage_cut;
    dworker.ruIndex = ruIndex;
    dworker.max_hang = max_hang;
    dworker.min_ovlp = min_ovlp;

    // (To guarantee that set_R_to_U's functionality can be
    // in threads (reallocation avoided beforehand) )
    if (ruIndex->len<n_read){
        fprintf(stderr, "[W::%s] Expanding R_to_U buffer, this might be unexpected go check code. old size %d, n_read %d\n", 
                __func__, (int)ruIndex->len, (int)n_read);
        ruIndex->index = (uint32_t*)realloc(ruIndex->index, n_read*sizeof(uint32_t));
        memset(ruIndex->index + ruIndex->len, 
                -1, 
                sizeof(uint32_t)*(n_read - ruIndex->len));
        ruIndex->len = n_read;
    }

    // mark things
    kt_for(asm_opt.thread_num, hit_contained_advance_callback, 
            &dworker, n_read);
    int tot_marked = 0;
    for (int i=0; i<n_read; i++){
        tot_marked += (int)todel2[i];
    }
    fprintf(stderr, "[T::%s] step1: used %.2f s, marked %d in bit array\n", 
            __func__, Get_T()-startTime, tot_marked);
    startTime = Get_T();

    // delete things
    // No need to search in the opposite direction when deleting, we have all the 
    // marks, so one pass for paf and one pass for coverage_cut should be enough.
    uint64_t tot_del = 0, tot_tried_reads=0;
    for (uint32_t i=0; i<n_read; i++){
        if (coverage_cut[i].del) continue;
        tot_tried_reads++;

        ma_hit_t_alloc *x = &sources[i];
        int just_del_all = todel2[i];
        if (just_del_all){
            for (int j=0; j<x->length; j++){
                if (!x->buffer[j].del) tot_del++;
                x->buffer[j].del = 1;
            }
        }else{
            // self is not deleted, but need to 
            // check whether any target is a dangling hit
            for (int j=0; j<x->length; j++){
                uint32_t tn = Get_tn(x->buffer[j]);
                if (todel2[tn]){
                    if (!x->buffer[j].del) tot_del++;
                    x->buffer[j].del = 1;
                }
            }
        }
    }
    free(todel2);
    fprintf(stderr, "[T::%s] step2: used %.2f s; deleted %" PRIu64 " (tried %" PRIu64 "reads))\n", 
            __func__, Get_T()-startTime, tot_del, tot_tried_reads);

    // update index
    startTime = Get_T();
    transfor_R_to_U(ruIndex);
    fprintf(stderr, "[T::%s] step2.5 used %.2f s\n", __func__, Get_T()-startTime); 

    // mark reads that have no non-contained neighbors and drop them
    startTime = Get_T();
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
            //fprintf(stderr, "[dbg::%s] del %d\n", __func__, (int)i);
            ret++;
            coverage_cut[i].del = 1;
        }
    }
    fprintf(stderr, "[T::%s] step3 used %.2f s\n", __func__, 
                Get_T()-startTime);

    fprintf(stderr, "[M::%s] dropped %d reads, used total of %0.2f s\n\n", __func__, ret, Get_T()-startTime0);
    return ret;
}
uint16_t get_i_from_paf_index(paf_ht_t *h, uint32_t qn, uint32_t tn, int which){
    // find the i in (reverse_)sources[tn].buffer[i] such that tn of it equals qn.
    // which: 0 for cis, 1 for trans
    // return (uint16_t)-1 if not found.
    khint_t key = paf_ht_get(h, ((uint64_t)qn)<<32 | tn );
    if (key!=kh_end(h)){
        uint16_t ret = kh_val(h, key);
        if ((ret&1) == which)
            return ret>>1;
    }
    return (uint16_t)-1;
}


paf_ht_t* hamt_index_pafs(ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,
                    ma_sub_t *coverage_cut,
                    long long n_reads){
    // Need to have quick access to indcies in pafs. Searching is too costly.
    // This produces one hashtable, pairing entries sources and reverse sources.
    // Key is packed read IDs: qn<<32|tn
    // value is: (index on the *current side*)<<1 | is_trans

    double T = Get_T();
    paf_ht_t *h = paf_ht_init();
    khint_t key;
    int absent;
    uint64_t hash, tmp;

    ma_hit_t_alloc *paf;
    ma_hit_t *v;
    uint32_t tn;
    uint64_t tot = 0;
    for (uint16_t is_trans=0; is_trans<=1; is_trans++){
        if (!is_trans) paf = sources;
        else paf = reverse_sources;
        if (!paf) continue;
        
        for (uint32_t qn=0; qn<n_reads; qn++){
            //if (coverage_cut[qn].del) continue;  // TODO: check this
            for (uint16_t i=0; i<paf[qn].length; i++){
                if (paf[qn].buffer[i].del) continue;
                v = &paf[qn].buffer[i];
                tn = Get_tn(*v);
                hash = ((uint64_t)qn)<<32 | tn; 
                key = paf_ht_put(h, hash, &absent);
                if (!absent){
                    fprintf(stderr, "[W::%s] overlap qn-tn pair "
                            "not unique? qn %d (sancheck %d) tn %d i %d. "
                            "hash upper 32 %" PRIu32 " lower 32 %" PRIu32 ".Dump: \n", 
                            __func__, (int)qn, (int)Get_qn(*v), (int)tn, (int)i, (uint32_t)(hash>>32), (uint32_t) hash);
                    for (int j=0; j<paf[qn].length; j++){
                        fprintf(stderr, "qn->tn: tn = %d\n", (int)Get_tn(paf[qn].buffer[j]));
                    }
                    for (int j=0; j<paf[tn].length; j++){
                        fprintf(stderr, "tn->qn: tn = %d\n", (int)Get_tn(paf[tn].buffer[j]));
                    }
                    continue;
                }
                kh_val(h, key) = i<<1 | is_trans;
                tot++;
            }
        }
    }
    fprintf(stderr, "[T::%s] inserted %" PRIu64 ", used %.2f s\n", 
            __func__, tot, Get_T()-T);
    return h;
}


struct paf_ct_v{
    int prefix_l;
    paf_ct_t **hs;
};
paf_ct_v *init_paf_ct_v(int prefix_l){
    paf_ct_v *H = (paf_ct_v*)calloc(1, sizeof(paf_ct_v));
    H->prefix_l = prefix_l;
    H->hs = (paf_ct_t**)calloc(1<<H->prefix_l, sizeof(paf_ct_t*));
    for (int i=0; i<1<<H->prefix_l; i++){
        H->hs[i] = paf_ct_init();
    }
    return H;
}
void destroy_paf_ct_v(paf_ct_v *H){
    for (int i=0; i<1<<H->prefix_l; i++){
        paf_ct_destroy(H->hs[i]);
    }
    free(H->hs);
    free(H);
}

typedef struct{
    int n_threads;
    int prefix_l;
    uint64_t prefix_mask;
    int sancheck_countermax;
    ma_hit_t_alloc *paf;
    uint64_t is_trans;
    uint32_t batchsize, start, n_reads;
    paf_ct_v *H;
}hamt_paf_ht_v_pl_t;  // pipeline struct
 
typedef struct{
    hamt_paf_ht_v_pl_t *pl;
    vu64_t *ahashes;  // (1<<prefix_l) arrays, not piggyback because not enough bits
    vu32_t *aindices; // (1<<prefix_l) arrays
}hamt_paf_ht_v_ps_t; // pipeline step struct

void destroy_hamt_paf_ht_v_ps_t(hamt_paf_ht_v_ps_t* s){
    for (int i=0; i<=s->pl->prefix_mask; i++){
        kv_destroy(s->ahashes[i]);
        kv_destroy(s->aindices[i]);
    }
    free(s->ahashes);
    free(s->aindices);
    free(s);
}
 
static void hamt_index_pafs_multithread_pipeline_worker(void *data, long prefixID, int tID){
    hamt_paf_ht_v_ps_t *s = (hamt_paf_ht_v_ps_t*)data;
    vu64_t *v = &s->ahashes[prefixID];
    vu32_t *vi = &s->aindices[prefixID];
    paf_ct_t *h = s->pl->H->hs[prefixID];
    khint_t key;
    int absent, prefix_l = s->pl->prefix_l;
    uint64_t sanhash;

    for (size_t i=0; i<v->n; i++){
        key = paf_ct_put(h, v->a[i], &absent);
        if (!absent){
            fprintf(stderr, "[E::%s] slot not empty, query is %" PRIu64 "\n", 
                    __func__, v->a[i]);
        }
        //kh_key(h, key) = v->a[i] ;
        kh_val(h, key) = vi->a[i];
    }
}
static void *hamt_index_pafs_multithread_pipeline(void *data, int step, void *in){
    hamt_paf_ht_v_pl_t *p = (hamt_paf_ht_v_pl_t*)data;
    double T0;
    int verbose = 0;
    if (step==0){
        uint32_t tot = 0, qn, tn;
        uint64_t hash, prefix;
        ma_hit_t_alloc *h;
        ma_hit_t *hh;
        hamt_paf_ht_v_ps_t *s = (hamt_paf_ht_v_ps_t*)calloc(1, sizeof(hamt_paf_ht_v_ps_t));
        s->pl = p;
        s->ahashes = (vu64_t*)calloc(1<<p->prefix_l, sizeof(vu64_t));
        s->aindices = (vu32_t*)calloc(1<<p->prefix_l, sizeof(vu32_t));
        for (int i=0; i<1<<p->prefix_l; i++){
            kv_init(s->ahashes[i]);
            kv_resize(uint64_t, s->ahashes[i], 128);
            kv_init(s->aindices[i]);
            kv_resize(uint32_t, s->aindices[i], 128);
        }
        // collect values
        uint32_t i;
        for (i=p->start; i<p->n_reads && i<p->start+p->batchsize; i++){
            h = &p->paf[i];
            int tmpl = h->length;
            if (h->length > UINT32_MAX/*p->sancheck_countermax*/){
                fprintf(stderr, "[W::%s] Too many targets at read#%d (%.*s): %" PRIu64 ". "
                        "Indexing will truncate wrt current sorting order.\n", 
                        __func__, (int)i, (int)Get_NAME_LENGTH(R_INF, i), 
                        Get_NAME(R_INF, i), (uint64_t)h->length
                        );
                tmpl = p->sancheck_countermax;
            }
            for (uint32_t j=0; j<tmpl; j++){
                hh = &h->buffer[j];
                //if (hh->del) continue; // note: commented out because ha(v0.13)'s function doesn't check this.

                qn = Get_qn(*hh);
                tn = Get_tn(*hh);
                hash = ((uint64_t)qn)<<32 | tn;
                prefix = hash & p->prefix_mask;
                //hash = (hash>>p->prefix_l)<<p->prefix_l | j;
                kv_push(uint64_t, s->ahashes[prefix], hash);
                kv_push(uint32_t, s->aindices[prefix], j);
                tot++;
            }
        }
        p->start = i;
        if (tot==0){
            destroy_hamt_paf_ht_v_ps_t(s);
        }else{
            if (verbose) fprintf(stderr, "[dbg::%s] step1 start=%d\n", __func__, (int)p->start);
            return s;
        }
    }else if (step==1){
        hamt_paf_ht_v_ps_t *s = (hamt_paf_ht_v_ps_t*)in;
        kt_for(s->pl->n_threads, hamt_index_pafs_multithread_pipeline_worker, s, 1<<s->pl->prefix_l);
        destroy_hamt_paf_ht_v_ps_t(s);
    }else{
        fprintf(stderr, "[E::%s] pipeline doesn't have this step, check code.\n", __func__);
        exit(1);
    }
    return 0;
}
paf_ct_v* hamt_index_pafs_multithread(ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,
                    ma_sub_t *coverage_cut,
                    long long n_reads, int batchsize, int n_threads ){
    double T = Get_T();

    ma_hit_t_alloc *paf;
    ma_hit_t *v;

    paf_ct_v *H = init_paf_ct_v(PAF_CT_PRE);

    hamt_paf_ht_v_pl_t pl;
    pl.n_threads = n_threads;
    pl.prefix_l = PAF_CT_PRE;
    pl.prefix_mask = (1<<PAF_CT_PRE)-1;
    //pl.sancheck_countermax = (1<<(PAF_CT_PRE-1)) -1;  // need 1bit for is_trans
    pl.sancheck_countermax = PAF_CT_PRE_M;
    pl.batchsize = batchsize;
    pl.n_reads= n_reads;
    pl.H = H;

    for (uint64_t is_trans=0; is_trans<=1; is_trans++){
        if (!is_trans) paf = sources;
        else paf = reverse_sources;
        if (!paf) continue;

        pl.paf = paf;
        pl.start = 0;
        pl.is_trans = is_trans;

        kt_pipeline(2, hamt_index_pafs_multithread_pipeline, &pl, 2);
    }

    return H;
}

int get_paf_index_ct(paf_ct_v *h, uint32_t qn, uint32_t tn){
    uint64_t hash = ((uint64_t)qn)<<32 | tn;
    paf_ct_t *hh = h->hs[hash&PAF_CT_PRE_M];
    khint_t key = paf_ct_get(hh, hash);
    if (key==kh_end(hh)){
        return -1;
    }else{
        //uint64_t tmp = PAF_CT_PRE_M & kh_key(hh, key);
        uint32_t tmp = kh_val(hh, key);
        return (int)tmp;
    }
}
int insert_paf_index_ct(paf_ct_v *h, uint32_t qn, uint32_t tn, uint32_t idx){
    int absent;
    int silent = 0;

    uint64_t hash = ((uint64_t)qn)<<32 | tn;
    paf_ct_t *hh = h->hs[hash&PAF_CT_PRE_M];
    khint_t key = paf_ct_put(hh, hash, &absent);
    if (!absent && !silent){
        fprintf(stderr, "[W::%s] inserted to non-empty slot. qn %d tn %d\n", 
                __func__, (int)qn, (int)tn);
    }
    //kh_key(hh, key) = ((hash)>>PAF_CT_PRE)<<PAF_CT_PRE | idx;
    kh_val(hh, key) = idx;
    return absent;
}



typedef struct {
    ma_hit_t_alloc *sources;
    const ma_hit_t_alloc* reverse_sources; 
    //paf_ht_t *paf_h;
    paf_ct_v *paf_h_cis;
    paf_ct_v *paf_h_trans;
    double *dTs_trans;  // per-thread timer 
    double *dTs_any;
}hamt_cleanweakhit_t;
static void hamt_clean_weak_ma_hit_t_worker(void *data, long jobID, int tID){  // Callback for kt_for
    hamt_cleanweakhit_t *s = (hamt_cleanweakhit_t*) data;
    uint32_t rID = jobID, qn, tn_weak, tn_strong, qs, qe;  
    ma_hit_t_alloc *h = &s->sources[rID];
    const ma_hit_t_alloc *hr;

    double dT_trans = 0, dT_any = Get_T();
    khint_t key;
    uint64_t hash;
    int absent;
    for (uint32_t i=0; i<h->length; i++){
        if (h->buffer[i].del) continue;
        if (h->buffer[i].ml==0){  // is weak overlap
            qn = Get_qn(h->buffer[i]);  // qn should be just rID
            tn_weak = Get_tn(h->buffer[i]);
            qs = Get_qs(h->buffer[i]);
            qe = Get_qe(h->buffer[i]);
            
            // (`inline int check_weak_ma_hit`. 
            // The following is need to be checked multiple times because we want to make
            // sure that the strong overlap supporting the deletion should  
            // explicitly come from implied other haplotype. This is doen by checking 
            // if tn_strong -> tn_weak is a trans overlap. It is not required 
            // that tn_strong -> tn_weak exists in the cis overlaps.)
            //
            // The check can't be reduced to a pre-built array, either do linear/binary search
            // or make a hashtable.
            int can_delete = 0, is_trans, idx;
            int san_jj1=-1, san_jj2=-1, tn1=-1, tn2=-1;

            for (int j=0; j<h->length; j++){
                if (j==i) continue;
                if (!h->buffer[j].del && 
                    h->buffer[j].ml==1 &&  // is strong
                    Get_qs(h->buffer[j]) <=qs && Get_qe(h->buffer[j])>=qe
                   ){ 
                    tn_strong = Get_tn(h->buffer[j]);
                    int is_found = 0;
                    double Ttmp = Get_T();

//                    // method1: search
//                    hr = &s->reverse_sources[tn_strong];
//                    for (int jj=0; jj<hr->length; jj++){
//                        if (hr->buffer[jj].del) continue;
//                        if (Get_tn(hr->buffer[jj])==tn_weak){
//                            is_found = 1;
//                            san_jj1 = jj;
//                            tn1 = tn_strong;
//                            break;
//                        }
//                    } 
                    // method2: hashtable
                    //hash = ((uint64_t)tn_strong)<<32 | tn_weak;
                    //key = paf_ht_get(s->paf_h, hash);
                    idx = get_paf_index_ct(s->paf_h_trans, tn_strong, tn_weak);
                    //if (key!=kh_end(s->paf_h) && (kh_val(s->paf_h, key)&1) ){
                    if (idx!=-1){
                        is_found = 1;
                        //san_jj2 = (int)(kh_val(s->paf_h, key) >>1);
                        san_jj2 = idx;
                        tn2 = tn_strong;
                    }

                    dT_trans+=Get_T()-Ttmp;
                    if (is_found){
                        can_delete = 1;
                        break;
                    }
                }
            }
//            if (san_jj1!=san_jj2){
//                fprintf(stderr, "[dbg::%s] jj1 %d jj2 %d, qn %d tn_weak %d, "
//                        "tn_strong1 %d tn_strong2 %d\n", 
//                        __func__, san_jj1, san_jj2, qn, tn_weak, tn1, tn2);
//            }

            if (!can_delete) continue;

            // now mark self overlap and the other way around
            int sancheck = 0;
//            hr = &s->sources[tn_weak];
//            for (int j=0; j<hr->length; j++){
//                if (Get_tn(hr->buffer[j])==qn) {
//                    hr->buffer[j].bl = 0;
//                    sancheck = 1;
//                    break;
//                }
//            }
            hash = ((uint64_t)tn_weak)<<32 | qn;
            //key = paf_ht_get(s->paf_h, hash);
            idx = get_paf_index_ct(s->paf_h_cis, tn_weak, qn);
            //if (key!=kh_end(s->paf_h)){
            if (idx!=-1){
                //uint16_t jj = kh_val(s->paf_h, key);
                uint16_t jj = (uint16_t)idx;
                //if (!(jj&1)){
                if (1/*!is_trans*/){
                    //s->sources[tn_weak].buffer[jj>>1].bl = 0;  // last bit indicates cis/trans
                    s->sources[tn_weak].buffer[idx].bl = 0;
                    sancheck = 1;
                }
            }
            if (!sancheck){
                fprintf(stderr, "[E::%s] qn->tn exists but tn->qn not found or is trans?"
                        "qn=%d tn=%d\n", 
                        __func__, (int)qn, (int)tn_weak);
                //exit(1);
            }else{  // the other way around was fine, go delete self
                h->buffer[i].bl = 0;
            }
        }

    }
    s->dTs_trans[tID] += dT_trans;
    s->dTs_any[tID] += Get_T()-dT_any;
}

/**
 * @func For each cis-overlap of each read, 
 *        if it has at least one strong overlap and the target read of 
 *        that strong overlap does not trans-overlap with the current target, 
 *        mark-delete the current overlap.
*/
int hamt_clean_weak_ma_hit_t2(ma_hit_t_alloc* const sources, 
                             ma_hit_t_alloc* const reverse_sources, 
                             ma_sub_t *coverage_cut, 
                             const long long n_reads,
                             paf_ct_v *paf_h_cis, 
                             paf_ct_v *paf_h_trans){
    // threaded `clean_weak_ma_hit_t`
    
    int verbose = 0;
    int n_treated;
    int n_threads = asm_opt.thread_num;
    uint32_t i, j, qn, tn;

    double T = Get_T();
    double T0 = T;

    // collect indexing
    //paf_ct_v *paf_h_cis = hamt_index_pafs_multithread(sources, 0, 
    //        coverage_cut, n_reads, 2048, asm_opt.thread_num);
    //paf_ct_v *paf_h_trans = hamt_index_pafs_multithread(0, reverse_sources, 
    //        coverage_cut, n_reads, 2048, asm_opt.thread_num);
    hamt_cleanweakhit_t data;
    data.sources = sources;
    data.reverse_sources = reverse_sources;
    data.dTs_trans = (double*)calloc(n_threads, sizeof(double));
    data.dTs_any = (double*)calloc(n_threads, sizeof(double));
    data.paf_h_cis = paf_h_cis;
    data.paf_h_trans = paf_h_trans;
    kt_for(n_threads, hamt_clean_weak_ma_hit_t_worker, &data, n_reads);

    double dT_trans=0, dT_any = 0;
    for (int i=0; i<n_threads; i++){
        dT_trans+=data.dTs_trans[i];
        dT_any+=data.dTs_any[i];
    }
    free(data.dTs_trans);
    free(data.dTs_any);
    fprintf(stderr, "[T::%s] step 1 used %.2f s. (cpu time %f s, trans check used %f s)\n", 
            __func__, Get_T()-T, dT_any, dT_trans);

    T = Get_T();
    // i think if we don't need to report n_treated, and that ha functions before 
    // this function call do not have unrealized changes of .bl makrings, 
    // the following block can be put inside the kt_for callback. 
    n_treated = 0;
    for (int i=0; i<n_reads; i++){
        for (j=0; j<sources[i].length; j++){
            if (sources[i].buffer[j].del) continue;
            if (sources[i].buffer[j].bl==0){
                sources[i].buffer[j].del = 1;
                n_treated++;
            }else{
                sources[i].buffer[j].del = 0;
            }
        }
    }
    //destroy_paf_ct_v(paf_h);
    fprintf(stderr, "[T::%s] step2 used %.2f s (should be short, otherwise read comments.\n", 
            __func__, Get_T()-T);

    fprintf(stderr, "[M::%s] treated %d, used %.2f s\n", __func__, 
                n_treated, Get_T()-T0);
    return n_treated;
}


typedef struct {
    ma_hit_t_alloc *paf;
    paf_ct_v *H;
    uint64_t *cnt1, *cnt2;  // for sancheck
}normalize_paf_s_t;

static void hamt_normalize_paf_worker(void *data, long jobID, int tID){
    normalize_paf_s_t *s = (normalize_paf_s_t*)data;
    uint32_t rID = jobID, qn, tn;
    uint32_t ql, ql_rev;
    ma_hit_t_alloc *h = &s->paf[rID];
    ma_hit_t *hh, *hh_rev;
    int is_del;

    uint64_t cnt1=0, cnt2 = 0;
    for (uint32_t i=0; i<h->length; i++){
        // (((ha doesn't require an unset .del of self here)))
        if (h->buffer[i].del) continue;

        hh = &h->buffer[i];
        qn = Get_qn(*hh);
        tn = Get_tn(*hh);

        // does the opposite exist?
        int j=get_paf_index_ct(s->H, tn, qn);
        if (j<0) {
            fprintf(stderr, "[E::%s] ovlp from the other way around not found."
                    "Will be fatal for later functions, check code now\n", __func__);
            exit(1);
        }

        is_del = 0;
        hh = &h->buffer[i];
        hh_rev = &s->paf[tn].buffer[j];
        if (hh->del || hh_rev->del) is_del = 1;

        ql =     Get_qe(*hh)     - Get_qs(*hh);
        ql_rev = Get_qe(*hh_rev) - Get_qs(*hh_rev);

        // always only update when qn<tn to avoid data race
        if (qn<tn){
            if (ql>=ql_rev){  // overwrite tn's
                set_reverse_overlap(hh_rev, hh);
                cnt1++;
            }else{  // overwrite self(shorter)
                set_reverse_overlap(hh, hh_rev);
                cnt2++;
            }
        }
        hh->del = is_del;
        hh_rev->del = is_del;
    }
    s->cnt1[tID]+=cnt1;
    s->cnt2[tID]+=cnt2;
}



typedef struct {
    ma_hit_t_alloc *paf;
    paf_ct_v *H;
    vu64_t *qni;
    uint64_t *cnt;
}symmetrize_paf_s_t;

static void hamt_symmetrize_paf_worker(void *data, long rID, int tID){
    symmetrize_paf_s_t *s = (symmetrize_paf_s_t*)data;
    ma_hit_t_alloc *h =  &s->paf[rID];
    uint32_t qn=rID, tn;
    int idx;
    uint64_t tmp;
    for (uint16_t i=0; i<h->length; i++){
        tn = Get_tn(h->buffer[i]);
        idx = get_paf_index_ct(s->H, tn, qn);
        if (idx<0){
            tmp = ((uint64_t)qn)<<32 | i;
            kv_push(uint64_t, s->qni[tID], tmp);
            s->cnt[tID]++;
        }
    }
}
static void hamt_symmetrize_by_del_paf_worker(void *data, long rID, int tID){
    symmetrize_paf_s_t *s = (symmetrize_paf_s_t*)data;
    ma_hit_t_alloc *h =  &s->paf[rID];
    uint32_t qn=rID, tn;
    int idx;
    uint64_t tmp;
    for (uint16_t i=0; i<h->length; i++){
        tn = Get_tn(h->buffer[i]);
        idx = get_paf_index_ct(s->H, tn, qn);
        if (idx<0){
            h->buffer[i].del = 1;
            s->cnt[tID]++;
        }
    }
}
void hamt_symmetrize_paf_add_new_entry(ma_hit_t_alloc *paf, paf_ct_v *H, 
                                        uint32_t qn, uint32_t qi){
    ma_hit_t ele;
    set_reverse_overlap(&ele, &paf[qn].buffer[qi]);
    paf[qn].buffer[qi].del = 1;
    ele.del = 1;
    uint32_t tn = Get_qn(ele);
    add_ma_hit_t_alloc(&paf[tn], &ele);
    insert_paf_index_ct(H, tn, qn, paf[tn].length-1);
}
void hamt_symmetrize_paf(ma_hit_t_alloc *paf, paf_ct_v *H, long long n_reads, int n_threads){
    // Symmetrize uses weirdly lot of memory? Switching to mark single entries
    // as deleted. (Was insert the counterpart then mark boths as deleted.)
    double T0 = Get_T();
    symmetrize_paf_s_t s;
    s.paf = paf;
    s.H = H;
    s.qni = (vu64_t*)calloc(n_threads, sizeof(vu64_t));
    s.cnt = (uint64_t*)calloc(n_threads, sizeof(uint64_t));
    for (int i=0; i<n_threads; i++){
        kv_init(s.qni[i]);
        kv_resize(uint64_t, s.qni[i], 64);
    }

    // collect asym entries
    kt_for(n_threads, hamt_symmetrize_by_del_paf_worker, &s, n_reads);
    for (int i=0; i<n_threads; i++){s.cnt[0]+=s.cnt[i];}
    fprintf(stderr, "[M::%s] step1, total %" PRIu64 " . used %.2f s\n",
            __func__, s.cnt[0], Get_T()-T0);

    //T0 = Get_T();
    //for (int i=0; i<n_threads; i++){
    //    radix_sort_ovhamt64(s.qni[i].a, s.qni[i].a+s.qni[i].n);
    //}
    //fprintf(stderr, "[M::%s] step1b, sorting used %.2f s (for cache efficiency, "
    //        "maybe doesn't matter. remove if takes a long time)\n", 
    //        __func__, Get_T()-T0);

    //// create new entries and update hashtables
    //T0 = Get_T();
    //uint32_t qn, qi;
    //for (int i=0; i<n_threads; i++){
    //    for (int j=0; j<s.qni[i].n; j++){
    //        qn = s.qni[i].a[j]>>32;
    //        qi = (uint32_t)s.qni[i].a[j];
    //        hamt_symmetrize_paf_add_new_entry(paf, H, qn, qi);
    //    }
    //}
    //fprintf(stderr, "[M::%s] step2, used %.2f s\n", __func__, Get_T()-T0);


    for (int i=0; i<n_threads; i++){
        kv_destroy(s.qni[i]);
    }
    free(s.qni);
    free(s.cnt);
}

void hamt_normalize_paf(ma_hit_t_alloc *paf, paf_ct_v *H, long long n_reads, int n_threads){
    // Rewrite hamt_normalize_ma_hit_t_single_side_advance
    // Use paf indexing.
    // Different from ha's function: the lack of counterpart ovlp should be
    //  already accounted for before calling this function. The ha's function
    //  does that and the normalization at the same time. This might lead to 
    //  minor differences? (not sure)

    double T0 = Get_T();    
    int verbose = 0;

    normalize_paf_s_t s;
    s.paf = paf;
    s.H = H;
    s.cnt1 = (uint64_t*)calloc(n_threads, sizeof(uint64_t));
    s.cnt2 = (uint64_t*)calloc(n_threads, sizeof(uint64_t));

    kt_for(n_threads, hamt_normalize_paf_worker, &s, n_reads);
    for (int i=1; i<n_threads; i++){
        s.cnt1[0] += s.cnt1[i];
        s.cnt2[0] += s.cnt2[i];
    }
    fprintf(stderr, "[M::%s] type1 %" PRIu64 ", type2 %" PRIu64 ", used %.2f s\n", 
            __func__, s.cnt1[0], s.cnt2[0], Get_T()-T0);
    free(s.cnt1);
    free(s.cnt2);
}

