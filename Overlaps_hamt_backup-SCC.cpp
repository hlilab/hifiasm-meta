// backup note:
//     Oct 14 2020
//     3 attemps of SCC, all not working correctly
//     i need SCC on bidirectional graph, not tricks.

#include <stdint.h>
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>
#include <assert.h>
#include "Overlaps.h"
#include "Process_Read.h"
#include "CommandLines.h" 
#include "Overlaps_hamt.h"
#include "ksort.h"

KDQ_INIT(uint64_t)
KDQ_INIT(uint32_t)
KRADIX_SORT_INIT(ovhamt64, uint64_t, uint64_t, 8)

//////////////////////////////////////////////////////////////////////
//                        debug functions                           //
//////////////////////////////////////////////////////////////////////
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
    for (int i=0; i<ug->u.n; i++){
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
    uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
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
    free(primary_flag);    
}

// void ma_usg_print_lite(const ma_ug_t *ug, asg_t *usg, All_reads *RNF, asg_t* read_g,
//                        int print_seq, const char* prefix, FILE *fp){

//     uint32_t *v2vu = (uint32_t*)calloc(usg->n_seq, sizeof(uint32_t));
//     uint8_t* primary_flag = (uint8_t*)calloc(read_g->n_seq, sizeof(uint8_t));
// 	uint32_t i, j, l, vu, v;

//     for (vu=0; vu<ug->u.n; vu++){
//         v2vu[ug->u.a[vu].start>>1] = vu;
//         v2vu[ug->u.a[vu].end>>1] = vu;
//     }

// 	char name[200];
// 	for (i = 0; i < usg->n_seq*2; ++i) { // the Segment lines in GFA
//         if (!usg->whitelist[i]){
//             continue;
//         }
//         vu = v2vu[i>>1];
//         sprintf(name, "%s%.6d(%.*s)", prefix, vu+1, (int)Get_NAME_LENGTH((*RNF), i>>1), Get_NAME((*RNF), i>>1));
// 		fprintf(fp, "S\t%s\t*\tLN:i:%d\tdp:f:1\n", name, 9);
        
//         fprintf(fp, "A\t%s\t%d\t%c\t%.*s\t%d\t%d\tid:i:%d\tHG:A:%c\n", name, l, "+",
//                      (int)Get_NAME_LENGTH((*RNF), i>>1), Get_NAME((*RNF), i>>1), 
//                      0, 10000, i>>1, "a");
// 	}
//     // for(i=0; i<usg->n_arc; i++){
//     //     uint32_t u = usg->arc[i].ul>>32, v = usg->arc[i].v;
// 	// 	fprintf(fp, "L\t%s%.6d(%.*s)\t%c\t%s%.6d%c\t%c\t%dM\tL1:i:%d\n", 
//     //                 prefix, (u>>1)+1, (int)Get_NAME_LENGTH((*RNF), i>>1), Get_NAME((*RNF), i>>1), "+",
//     //                 prefix,	(v>>1)+1, (int)Get_NAME_LENGTH((*RNF), i>>1), Get_NAME((*RNF), i>>1), "+", ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
// 	// }
//     free(primary_flag);   
//     free(v2vu);
// }


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
    return nv++;
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
    return nv++;
}

int hamt_asgarc_util_available_vtx_of_node(asg_t *g, uint32_t nodeID, uint32_t *v){
    // given a node ID, store either of the vertices that has at least 1 arc to v
    // return 1 if found, 0 otherwise
    if (asg_arc_n(g, nodeID<<1)>0){
        *v = nodeID<<1;
        return 1;
    }else if (asg_arc_n(g, nodeID<<1 | 1)>0){
        *v = nodeID<<1 | 1;
        return 1;
    }else{
        return 0;
    }
}


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

        for (i=0; i<nv; i++)
            if (!av[i].del) 
                {v = av[0].v; 
                acc += (uint32_t) av[0].ul;
                break;}
        if (hamt_asgarc_util_countPre(sg, v, 0, 0)!=1) return 0;  // backward bifur 
    }
    return 1;
}


int hamt_asgutil_detect_dangling_circle(const asg_t *g, uint32_t v0, uint32_t w, uint8_t *buf_traversal){
    // EXPERIMENTAL and DEPRECATED; use strongly connect components routines instead!
    // (defaulted to ignore any arc/vertex that has been deleted)
    // detect if the given vertex is the start of a circle (search in one direction), even a very large one (*i dont want to abuse detect_single_path so here it is)
    // REQUIRE that w is not in the circle, reagardless of direction
    // probably expensive; only intend to be used with circle-aware coverage-based arc drop
    
    // return 0 if no, otherwise 1. Won't give the size of that circle, nor structure info.

    // ~~~~~ debug ~~~
    char str_buf[100];
    uint64_t str_len;
    //~~~~~~~~~~~~~~

    uint32_t v = v0, nv;
    asg_arc_t *av;

    // asg32_v queue = {0, 0, 0};
    // kv_push(uint32_t, queue, v0);
    kdq_t(uint32_t) *queue = kdq_init(uint32_t);
    queue->count = 0;
    kdq_push(uint32_t, queue, v0);
    // int early_termination = 0;

    // while (queue.n>0){
    while (kdq_size(queue) != 0){
        // fprintf(stdout, "queue: %d, ", (int)kdq_size(queue));fflush(stdout);
        v = *kdq_pop(uint32_t, queue);// v = kv_pop(queue);
        buf_traversal[v] = 1;  // need this because we don't know if there's circles in side of circles, i think
        av = asg_arc_a(g, v);
        nv = asg_arc_n(g, v);

        str_len = R_INF.name_index[(v>>1)+1]-R_INF.name_index[v>>1]; 
        memcpy(str_buf, R_INF.name+R_INF.name_index[v>>1], str_len);
        str_buf[str_len] = '\0';
        // fprintf(stdout, "at %s\n", str_buf);

        for (int i=0; i<(int)nv ;i++){
            str_len = R_INF.name_index[(av[i].v>>1)+1]-R_INF.name_index[av[i].v>>1]; 
            memcpy(str_buf, R_INF.name+R_INF.name_index[av[i].v>>1], str_len);
            str_buf[str_len] = '\0';

            if (av[i].del || g->seq[av[i].v>>1].del) {
                continue;
            }
            if (av[i].v>>1 == w>>1) {
                // w must not be in the loop
                // NOTE do not break or return 0 here! that's wrong.
                continue;
            }
            if (av[i].v == v0){  // NOTICE: requires that it's the same direction
                // fprintf(stdout, "  -> %s, HIT\n", str_buf);
                // free(queue.a);
                kdq_destroy(uint32_t, queue);
                return 1;
            } else{
                if (buf_traversal[av[i].v]==0){
                    // fprintf(stdout, "  -> %s, push\n", str_buf);
                }
            }
            kdq_push(uint32_t, queue, av[i].v);// kv_push(uint32_t, queue, av[i].v);
        }
        // early_termination++;
        // if (early_termination>10000){
        //     fprintf(stdout, "  -> took too long, terminate.\n");
        //     break;
        // }
        fflush(stdout);
    }
    fflush(stdout);
    // free(queue.a);
    kdq_destroy(uint32_t, queue);
    return 0;
    
}


int hamt_asgarc_util_get_the_one_target(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc){
    // v as exactly one target (sucessor), given the `include_del_seq` and `include_del_arc` (NOT tested in this function)
    // get that vertex's index in av (NOT vertex ID)
    int nv = asg_arc_n(g, v);
    if (nv==1) {
        return 0;
    }
    asg_arc_t *av = asg_arc_a(g, v);
    int i;
    for (i=0; i<(int)nv; i++){
        if (!include_del_arc && av[i].del) {continue;}
        if (!include_del_seq && g->seq[av[i].v>>1].del) {continue;}
        return i;
    }
    fprintf(stderr, "[E::%s] unexpected.\n", __func__); exit(1);
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
        wid = hamt_asgarc_util_get_the_one_target(g, v, include_del_seq, include_del_arc);
        av = asg_arc_a(g, v);
        w = av[wid].v;
        if (w>>1 == v0>>1){
            return -1;  // loop
        }

        nb_bwd =  hamt_asgarc_util_countPre(g, w, include_del_seq, include_del_arc);
        assert(nb_bwd!=0);
        if (nb_bwd>1) {break;}  // w has more than 1 predecessor

        v = w;
        if (max_arcs>0) {
            *traversal_dist = *traversal_dist +1;
            if (*traversal_dist>max_arcs) {return 0;}
        }else {
            *traversal_dist += (int) ((uint32_t)av[wid].ul);
            if (*traversal_dist>max_length) {return 0;}
        }

        if (a){
            kv_push(uint64_t, *a, (uint64_t)w);
        }
    }
    return 1;  // encountered non-single edge

}

int hamt_asgarc_util_checkSimpleBubble(asg_t *g, uint32_t v0, uint32_t *dest){
    // FUNC
    //     check if v0 is the start of a simple bubble (i.e. a 1-vertex bubble)
    // RETURN
    //     0 if no
    //     1 if yes
    if (hamt_asgarc_util_countSuc(g, v0, 0, 0)!=2){
        return 0;
    }
    uint32_t w[2], x, nv;  // x is the next-next vertex (if exists)
    asg_arc_t *av = asg_arc_a(g, v0);
    int idx = 0;
    
    for (uint32_t i=0; i<asg_arc_n(g, v0); i++){
        if (av[i].del) {continue;}
        if (g->seq[av[i].v>>1].del) {continue;}
        w[idx] = av[i].v;
        idx++;
    }
    assert(idx==2);
    // check if the next-next vertex is the end of a bubble
    if (hamt_asgarc_util_countSuc(g, w[0], 0, 0)!=1 || hamt_asgarc_util_countSuc(g, w[1], 0, 0)!=1){  // need to end in only one vertex
        return 0;
    }
    if (hamt_asgarc_util_get_the_one_target(g, w[0], 0, 0)!= hamt_asgarc_util_get_the_one_target(g, w[1], 0, 0)){  // and it's the same vertex
        return 0;
    }
    x = hamt_asgarc_util_get_the_one_target(g, w[0], 0, 0);
    if (hamt_asgarc_util_countPre(g, x, 0, 0)!=2){  // and the "end of the simple bubble" must not have other incoming sources
        return 0;
    }
    if (dest){
        *dest = x;
    }
    return 1;
}

int hamt_asgarc_util_countSinglePath_allowSimpleBubble(asg_t *g, uint32_t v0, int max_arcs){
    // REQUIRE
    //     graph is clean of loops!
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
            w = hamt_asgarc_util_get_the_one_target(g, v, 0, 0);
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
        }
        if (dist>max_arcs){
            return 0;
        }
    }
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
        wid = hamt_asgarc_util_get_the_one_target(g, w, 0, 0);
        w = asg_arc_a(g, w)[wid].v;
        if (hamt_asgarc_util_countPre(g, w, 0, 0)!=1){
            break;
        }
    }
    cov = cov/n;
    return (int) (cov+0.499);
}

int hamt_asgutil_calc_SCC_cov(asg_t *g, uint32_t v0, int *SCC_labels){
    // FUNC
    //     given v, find the kmer-based overage of reads of the same SCC group
    int label = SCC_labels[v0>>1];
    float cov = 0;
    int n = 0;
    for (uint32_t v=0; v<g->n_seq; v++){  // iterate over nodes, not vertices
        if (SCC_labels[v]==label){
            cov+=R_INF.median[v];
            n++;
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
    has_branching = hamt_asgarc_util_countSinglePath_allowSimpleBubble(g, v, 10000);  // 10k is the number of vertices, not bp
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
        // a->ul = ((uint64_t)a0->v)<<32 | (a0->ul>>32);
        a->ul = ((uint64_t)a0->v)<<32 | (uint64_t)((uint32_t)a0->ul);
        a->v = a0->ul>>32;
    }
    
    asg_arc_sort(g);
    g->is_srt = 1;
    asg_arc_index(g);

    return g;

}

void hamt_asggraph_util_DFS(asg_t *g0, uint64_t *ts_curr, uint64_t *ts_order, int *SCC_marking, int is_use_reverse){
    // REQUIRE
    //     assuming all .del marks have been applied
    // PAR
    //     `ts_curr` is an array that's allocated by the caller
    //        - ts's length is /*n_seq*/ total_reads (sparse but whatever)
    //        - for each entry, upper 32 bits is the start time, lower 32 bits is the end time
    //     `is_use_reverse`: if set to 1, reverse all arcs (i.e. use predecessor(s) instead of target(s) of an vertex)
    //                       requires *ts_order to be not NULL
    //                       (for SCC's 2nd traversal)
    //      `SCC_marking`: if not NULL and `is_use_reverse` is set, will store SCC labeling 
    // FUNC 
    //     do DFS on assembly graph g, seeding order is either based on vertexID order, or based on ts_order (SCC 2nd traversal)
    // SIDE EFFECT
    //     store traversal time to ts


    // original or transpose?
    asg_t *g;
    if (!is_use_reverse){
        g = g0;
    }else{
        g = hamt_asggraph_util_gen_transpose(g0);
    }

    // aux
    fprintf(stderr, "[D::%s] (entered).\n", __func__); fflush(stderr);
    double startTime = Get_T();
    int *color = (int*)calloc(R_INF.total_reads, sizeof(int));  // 0 for white, 1 for gray, 2 for black; -1 for orphan
    stacku32_t DFSstack;
    stacku32_init(&DFSstack);
    uint32_t nb_node_finished = 0;
    uint64_t DFStime = 0;
    uint32_t nv;
    asg_arc_t *av;
    int SCC_label = 1;
    uint32_t v, w;

    // mark orphan nodes
    // note (not sure):
    //   transitive reduction does not remove node/seq
    //   so here we need to test if a vertex has any arc left
    //   if not, we mark the node as -1 and don't check it in the 2nd DFS
    int ret_status;
    for (uint32_t i=0; i<R_INF.total_reads; i++){
        ret_status = hamt_asgarc_util_available_vtx_of_node(g, i, &v);
        if (!ret_status){  // node has no arc on either directions
            color[i] = -1;
            if (ts_curr){
                ts_curr[i] = (uint64_t) 0 | i<<1;
            }
            if (SCC_marking){
                SCC_marking[i] = 0;
            }
        }
    }

    // init
    uint32_t tmp_i;
    for (tmp_i=0; tmp_i<R_INF.total_reads; tmp_i++){
        if (!ts_order){
            if (color[tmp_i]==-1){continue;}
            assert(hamt_asgarc_util_available_vtx_of_node(g, tmp_i, &v));
            fprintf(stderr, "[debug::DFSinit] not given ordering, 1st vertex is %d (%d arcs).\n", (int)v, (int)asg_arc_n(g, v));
            break;
        }else{
            v = (uint32_t)ts_order[R_INF.total_reads-tmp_i-1];
            if (color[v>>1]==-1){continue;}
            fprintf(stderr, "[debug::DFSinit] given ordering, 1st vertex is %d (%d arcs).\n", (int)v, (int)asg_arc_n(g, v));
            break;
        }
    }
    if (tmp_i==R_INF.total_reads){
        fprintf(stderr, "ERROR, can't find the 1st vertex to start DFS.\n");
        exit(1);
    }
    color[v>>1] = 1;
    stacku32_push(&DFSstack, v);
    if (SCC_marking) {SCC_marking[v>>1] = SCC_label;}

    // init timestamp array (packing nodeID for radix sort)
    // if (ts_curr){
    //     uint64_t idx;
    //     // for (uint64_t i=0; i<g->n_arc; i++){  // is this correct? use n_seq?
    //     //     idx = g->arc[i].ul>>32;
    //     //     fprintf(stderr, "wut\t%d\n", (int)(idx>>1));
    //     //     ts_curr[idx>>1] = idx;
    //     // }
    //     for (uint64_t i=0; i<g->n_seq; i++){  // is this correct? use n_seq?
    //         ts_curr[i] = i;
    //     }
    // }

    // DFS
    //////////////////////////////////////////////////////////
    // not the recursive version,                           //
    // probably not the best way to check the color i guess?//
    // but please optimize later... not right now           //
    //////////////////////////////////////////////////////////
    while (1){
        DFStime++;
        // fprintf(stderr, ">>> DFStime %d, stack size is %d\n", (int)DFStime, DFSstack.n);
        if (stacku32_pop(&DFSstack, &v)){
            // fprintf(stderr, "    (node ID is %d)\n", (int)(v>>1));
            // fprintf(stderr, "    (node is %.*s)\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1));
            ///////////////////////////////
            //           debug           //
            // // check all colors
            // int sanity[4] = {0, 0, 0, 0};
            // for (uint32_t i=0; i<R_INF.total_reads; i++){
            //     sanity[color[i]+1]++;
            // }
            // fprintf(stderr, "     (orphan: %d, white: %d, gray: %d, black: %d\n)", sanity[0], sanity[1], sanity[2], sanity[3]);
            ///////////////////////////////


            // (if this node isn't finished, it'll be pushed right back in)

            // current node has been exhausted
            if (color[v>>1]==2){
                // fprintf(stderr, "continue, bcs node color\n");
                DFStime--;
                continue;
            }

            // examine the current node
            nv = asg_arc_n(g, v);
            av = asg_arc_a(g, v);
            // fprintf(stderr, "    (nv is %d)\n", (int)nv);
            if (nv==0){  // is a leaf node
                if (ts_curr) {ts_curr[v>>1] = DFStime<<32 | v;}
                color[v>>1] = 2;
                nb_node_finished++;
                // fprintf(stderr, "node is a leaf, flipping\n");
            }else{  // has child(ren)
                // check if all children are black
                // fprintf(stderr, "node is NOT a leaf\n");
                int is_exhausted = 1;
                for (uint32_t i=0; i<nv; i++){
                    if (color[av[i].v>>1]==0){
                        is_exhausted = 0;
                        break;
                    }
                }
                if (!is_exhausted){  // there's still unvisited child nodes
                    // fprintf(stderr, "node's children are not finished\n");
                    // put the current vertex back
                    stacku32_push(&DFSstack, v);
                    // push child nodes
                    for (uint32_t i=0; i<nv; i++){
                        w = av[i].v;
                        if (color[w>>1]==0){  // node hasn't been visited
                            // fprintf(stderr, "PUSHING one node\n");
                            color[w>>1] = 1;
                            stacku32_push(&DFSstack, w);
                            if (SCC_marking) {SCC_marking[w>>1] = SCC_label;}
                        }
                    }
                }else{  // subtree finished, update the current node's status
                    // fprintf(stderr, "node's children all finished, flipping\n");
                    if (ts_curr) {ts_curr[v>>1] = DFStime<<32 | v;}
                    color[v>>1] = 2;
                    nb_node_finished++;
                }
            }
        }else{  // stack is empty
            ///////////////////////////////
            //           debug           //
            // check all colors
            int sanity[4] = {0, 0, 0, 0};
            for (uint32_t i=0; i<R_INF.total_reads; i++){
                sanity[color[i]+1]++;
            }
            // assert(sanity[1]==0);  // there shouldn't be any gray node
            // assert(sanity[2]==nb_node_finished);  // the counter should reflect the number of black nodes             
            // fprintf(stderr, "DFS debug: stack emptied, nb_node_finished==%d\n", nb_node_finished);fflush(stderr);
            // fprintf(stderr, "           (orphan: %d, white: %d, gray: %d, black: %d)\n", sanity[0], sanity[1], sanity[2], sanity[3]);
            ///////////////////////////////

            if (nb_node_finished==R_INF.total_reads){  // check if there's any node left untouched
                // fprintf(stderr, "FINISHED, smooth.\n"); fflush(stderr);
                break;
            }else{  // get a new seed
                // fprintf(stderr, "~GET A NEW SEED.~\n"); fflush(stderr);
                SCC_label++;  // update SCC marking
                if (!ts_order){  // use the random order from vertex id
                    uint32_t i;
                    for (i=0; i<R_INF.total_reads; i++){
                        if (color[i]==0){
                            // v = i<<1 | 0;
                            assert(hamt_asgarc_util_available_vtx_of_node(g, i, &v));
                            break;
                        }
                    }
                    if (i==R_INF.total_reads){
                        // fprintf(stderr, "FINISHED, was searching but no node with arcs available left.\n");
                        break;
                    }
                }else{  // use the specified order
                    uint32_t j=0;
                    int i;
                    for ( i=R_INF.total_reads-1; i>=0; i--){
                        j = (uint32_t)ts_order[i];
                        if (color[j>>1]==0){
                            v = j;  // don't call hamt_asgarc_util_available_vtx_of_node! i think
                            break;
                        }
                    }
                    if (i==-1){
                        // fprintf(stderr, "FINISHED, was searching but no node with arcs available left.\n");
                        break;
                    }
                }
                color[v>>1] = 1;
                stacku32_push(&DFSstack, v);
                if (is_use_reverse && SCC_marking) {SCC_marking[v>>1] = SCC_label;}
            }
        }
    }
    stacku32_destroy(&DFSstack);
    free(color);
    if (is_use_reverse){
        asg_destroy(g);  // destroy the transposed graph
    }
    fprintf(stderr, "[D::%s] (exited).\n", __func__); fflush(stderr);
    fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
}

int *hamt_asggraph_util_SCC(asg_t *g){
    // REQUIRE
    //     assuming all .del marks have been applied
    // FUNC
    //     return labels of strongly connected components for nodes of *g

    fprintf(stderr, "[D::%s] (entered).\n", __func__); fflush(stderr);
    double startTime = Get_T();

    uint64_t *ts1 = (uint64_t*)calloc(R_INF.total_reads, sizeof(uint64_t));  // upper 32 bits: finish timestamp, lower 32 bits: node ID
    int *SCC_labels = (int*)calloc(R_INF.total_reads, sizeof(int));

    // 1st DFS
    hamt_asggraph_util_DFS(g, ts1, NULL, NULL, 0);

    // get the ordering: finish timestamp as of in 1st DFS
    radix_sort_ovhamt64(ts1, ts1+R_INF.total_reads);
    for (uint32_t i=0; i<R_INF.total_reads; i++){
        if ((ts1[i]>>32)==0){continue;}  // orphan node
        // fprintf(stderr, "radix\t%d\t%.*s\n", (int)(ts1[i]>>32), (int)Get_NAME_LENGTH(R_INF, (uint32_t)ts1[i]>>1), Get_NAME(R_INF, (uint32_t)ts1[i]>>1));
    }

    // 2nd DFS
    hamt_asggraph_util_DFS(g, NULL, ts1, SCC_labels, 1);


    free(ts1);
    if (VERBOSE>=1){
        fprintf(stderr, "[D::%s] (exited).\n", __func__); fflush(stderr);
        fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    }
    return SCC_labels;
}

void hamt_asgarc_SCCcovCut(asg_t *g){
    // mark strongly connected components, and cut the arcs with differential coverages that connect two SCCs

    asg_cleanup(g);
    asg_symm(g);
    fprintf(stderr, "[debug::%s], n_seq is %d, n_arc is %d, R_INF.total_reads is %d.\n", __func__, (int)g->n_seq, (int)g->n_arc, (int)R_INF.total_reads);

    // aux
    int n_reduced = 0;

    int *SCC_labels = hamt_asggraph_util_SCC(g);

    if (VERBOSE>=1){
        // debug
        int max_SCC_label = 1;
        for (int i=0; i<R_INF.total_reads; i++){
            if (SCC_labels[i]>max_SCC_label){
                max_SCC_label = SCC_labels[i];
            }
        }
        fprintf(stderr, "[D::%s] asg has %d strongly conneccted components.\n", __func__, max_SCC_label);
        write_debug_SCC_labeling(g, SCC_labels, asm_opt.output_file_name);
    }

    uint32_t v, w, n_vtx = g->n_seq*2;
    uint32_t nv, nv_;
    asg_arc_t *av, *av_;
    int handle_kmer_cov, SCC_kmer_cov;
    for (v=0; v<n_vtx; v++){
        nv = asg_arc_n(g, v);
        if (nv<1){  // note: allow the vertex in question to have more than 1 target
            continue;
        }
        nv_ = asg_arc_n(g, v^1);  // check backwards
        if ( nv_==0 || nv_>2 || (nv_==2 && !hamt_asgarc_util_checkSimpleBubble(g, v^1, NULL)) ){
            // no predecessor / multi incoming sources / 2 incoming sources but it's not a simple bubble
            continue;
        }
        av = asg_arc_a(g, v);

        for (int idx_w=0; idx_w<nv; idx_w++){
            w = av[idx_w].v;  // current read and the target read shall be in different SCCs
            if (SCC_labels[v>>1]==SCC_labels[w>>1]){
                continue;
            }
            // nv_ = asg_arc_n(g, w^1);  // check the target read backwards: if it's a circle, at least one predecessor of the target read should be in the same SCC with the target read
            // av_ = asg_arc_a(g, w^1);
            // int is_looking_like_a_circle = 0;
            // for (int idx_p=0; idx_p<nv_; idx_p++){  // check the target read's predecessor(s)
            //     if (SCC_labels[w>>1]==SCC_labels[av_[idx_p].v>>1]) {
            //         is_looking_like_a_circle = 1;
            //         break;
            //     }
            // }
            // if (!is_looking_like_a_circle){
            //     continue;
            // }
            
            // check coverage
            handle_kmer_cov = hamt_asgutil_calc_singlepath_cov(g, v);
            SCC_kmer_cov = hamt_asgutil_calc_SCC_cov(g, w, SCC_labels);
            fprintf(stderr, "[SCCdbg] handle %.*s | %d, targetSCC %d\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), handle_kmer_cov, SCC_kmer_cov);
            if ((float)SCC_kmer_cov / handle_kmer_cov <3){  // insufficient coverage diff
                continue;
            }else{
                if (SCC_kmer_cov<10){  // the SCC is not very well covered, don't touch
                    continue;
                }else{
                    // cut the arc
                    av[0].del = 1;
                    n_reduced++;
                }
            }
        }
    }

    if (n_reduced){
        asg_cleanup(g);
        asg_symm(g);
    }
    free(SCC_labels);

}


void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length){
    // rationale: tips could hinder simple bubble detection, and hamt doesn't need to be too careful with both of them (i think) (Sep 17)  
    double startTime = Get_T();
    uint32_t n_vtx = g->n_seq*2, v;
    uint32_t n_del;
    int pre, suc, l, ret;

    asg64_v a = {0,0,0};

    while (1){
        n_del = 0;
        // drop short tips
        for (v=0; v<n_vtx; v++){
            if (g->seq[v>>1].del) {continue;}
            pre = hamt_asgarc_util_countPre(g, v, 0, 0);
            suc = hamt_asgarc_util_countSuc(g, v, 0, 0);
            if (pre==0 && suc==1){  // tip of the other way around will be address later (or is already addressed) since we are iterating all vertices
                // traverse and see if single path down the way is short enough
                ret = hamt_asgarc_util_countSinglePath(g, v, 0, 0, max_arcs, max_length, &a, &l);
                // if (ret>=0){  // not a very long tip, also not a loop, within the search range
                if (ret>0){  // not a loop, and found the branching point within search range (otherwise it's either a really long tip, or we exhausted the whole utg)
                    ret = hamt_asgutil_is_tigSafeToErode(g, v);  // search again, this time also allow simple 1-vertex bubbles
                    if (ret){
                        // remove the whole tig
                        g->seq[v>>1].del = 1;
                        for (int i=0; i<(int)a.n; i++){
                            g->seq[(uint32_t)a.a[i]>>1].del = 1;
                            n_del+=1;
                        }
                    }
                }
            }
        }
        if (n_del>0){
            asg_cleanup(g);
            asg_symm(g);
        }

        // drop simple circles
        n_del += asg_arc_del_simple_circle_untig(NULL, NULL, g, 100, 0);

        // pop bubbles
        // (comment in ha was "remove isoloated single read")
        n_del += asg_arc_del_single_node_directly(g, asm_opt.max_short_tip, sources);  // includes graph cleanups

        // check
        if (n_del==0)
            {break;}
        else{
            if (VERBOSE>=1){
                fprintf(stderr, "[M::%s] removed %d locations\n", __func__, n_del);
                fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
            }
        }
    }
    free(a.a);

}

// DO NOT USE (this was before implementing asg SCC)
void hamt_asg_arc_del_by_readcov_circle_aware(asg_t *g){
    // drop an arc if 
    //    - the vertex in focus is not in a bubble
    //    - the vertex in focus is low coverage, 
    //    - the node only has only one arc in this direction
    //    - target and its targets on one side are ok/high coverage ones, AND they form a loop that does not involve the current read
    double startTime = Get_T();
    fprintf(stdout, "entering %s\n", __func__); fflush(stdout);
    int cnt = 0;
    uint32_t n_vtx = g->n_seq*2, v, v0;
    int cov, cov2;
    asg_arc_t *av;
    uint8_t *buf_traversal = (uint8_t*)calloc(n_vtx, sizeof(uint8_t));
    char str_buf[500]; uint64_t str_len; // debug
    for (v=0; v<n_vtx; v++){       
        str_len = R_INF.name_index[(v>>1)+1]-R_INF.name_index[v>>1]; 
        memcpy(str_buf, R_INF.name+R_INF.name_index[v>>1], str_len);
        str_buf[str_len] = '\0';
        // if (strcmp("m54337_181215_161949/63897767/ccs/worker3_136534", str_buf)!=0){continue;}

        if (hamt_asgarc_util_countSuc(g, v, 0, 0)!=1) {  // only check in the current direction; the other direction will be address since we iterate through all the vertices
            continue;
        }

        // get the "unitig" (the long single path in reversed direction)
        cov = hamt_asgutil_calc_singlepath_cov(g, v);
        fprintf(stdout, "m\t%s\t%d\n", str_buf, cov);
        // if (cov>30){
        //     fprintf(stdout, "median, %s\t%d\n", str_buf, cov);
        //     continue;
        // }

        fprintf(stdout, ">%s\n", str_buf);

        // (^~~~ maybe we can be more permissive, e.g. ignore if there's long tips connected to this vertex)
        v0 = asg_arc_a(g, v)[hamt_asgarc_util_get_the_one_target(g, v, 0, 0)].v;  // the single sucessor node
        if (hamt_asgarc_util_countPre(g, v0, 0, 0)<=1) {  // can't possibly be a circle
            fprintf(stdout, "=%s\tskip.defiNotCircle\n", str_buf); fflush(stdout);
            continue;
        }
        
        // if (R_INF.median[v0>>1] < (R_INF.median[v>>1]+1)*3) {  // require coverage difference
        cov2 = hamt_asgutil_calc_singlepath_cov(g, v0);
        if (cov*3 > cov2 || cov>50) {
            fprintf(stdout, "nodiff, %s\t%d\t%d\n", str_buf, cov, hamt_asgutil_calc_singlepath_cov(g, v0));
        }  
        memset(buf_traversal, 0, n_vtx*sizeof(uint8_t));
        if (hamt_asgutil_detect_dangling_circle(g, v0, v, buf_traversal)){ // check if it's in a (large) circle (allow any bubble/tangle/bifurcation/tip along the way)
            av = asg_arc_a(g, v);
            for (int i=0; i<(int)asg_arc_n(g, v); i++){  // (actually there's only 1 remaining arc but its index is not readily available, so.)
                av[i].del = 1;
                cnt++;
            }
            fprintf(stdout, "=%s\tYES\n", str_buf);fflush(stdout);
        }else{
            fprintf(stdout, "=%s\tnotcirc\n", str_buf);fflush(stdout);
        }
    }

    if (cnt>0) {
		asg_cleanup(g);
		asg_symm(g);
	}
    if (VERBOSE>=1){
        fprintf(stderr, "[M::%s] removed %d arcs\n", __func__, cnt);
        fprintf(stderr, "[M::%s] takes %0.2f s\n\n", __func__, Get_T()-startTime);
    }
    free(buf_traversal);

}




///////////////////////////////////////////////////////////////////////////////////////
//                  higher-level routines (unitig graph)                             //
///////////////////////////////////////////////////////////////////////////////////////

/*
                             NOTES
  - ug->g is a string graph that was built during unitig graph generation,
          i.e. it is not the "working" string graph, and it's not indexed by asg_arc_index(). (Oct 7 2020)

*/

void hamt_utg_redoug(ma_ug_t *ug, asg_t *sg){
	if (ug == 0) return;
    ma_ug_destroy(ug);
    ug = ma_ug_gen(sg);
}


ma_ug_t *hamt_utgutil_transpose(ma_ug_t *ug, asg_t *g){
    // NOTE: is only intended to be used in SCC routines.
    ma_ug_t *ret = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ret->g = hamt_asggraph_util_gen_transpose(g);  // note that this is the transpose of the "working" string graph, not the transpose of ug->g in ha's terms
    ret->u.m = ug->u.m;
    ret->u.n = ug->u.n;
    ret->u.a = (ma_utg_t*)calloc(ug->u.n, sizeof(ma_utg_t));

    ma_utg_t *h, *hh;
    for (int i=0; i<ret->u.n; i++){
        h = &ug->u.a[i];
        hh = &ret->u.a[i];
        hh->len = h->len; hh->circ = h->circ;
        // hh->start = h->end^1;
        // hh->end = h->start^1;
        hh->start = h->start;
        hh->end = h->end;


        // OMITING: h->a (list of reads), h->m, h->n
        // OMITING: h->s (unitig sequence)
    }
    return ret;
}


int test_utg_direction_require_target(ma_ug_t *ug, asg_t *sg, uint32_t vu, uint8_t *sanmask){
    // sanmask is an array that indicates whether a given vertex is a valid start/end vertex for the unitigs (vertex, not seg)
    // return 0 if default direction (start is start)
    //        1 if reverse (end is start)
    //        -1 if orphan
    int flag1=0, flag2 = 0;
    uint32_t nv1 = asg_arc_n(sg, ug->u.a[vu].start^1);
    uint32_t nv2 = asg_arc_n(sg, ug->u.a[vu].end^1);
    asg_arc_t *av1 = asg_arc_a(sg, ug->u.a[vu].start^1);
    asg_arc_t *av2 = asg_arc_a(sg, ug->u.a[vu].end^1);

    if (nv1>0){
        for (int i=0; i<nv1; i++){
            if (sanmask[av1[i].v]){
                flag1 = 1;
                break;
            }
        }
    }
    if (nv2>0){
        for (int i=0; i<nv2; i++){
            if (sanmask[av2[i].v]){
                flag2 = 1;
                break;
            }
        }
    }

    if (nv2>0 && flag2){  // start is start, end is end
        fprintf(stderr, "~test direction~ utg%.6d , 0\n", vu+1);
        return 0;
    }else if (nv1>0 && flag1){  // start is end, end is start
        fprintf(stderr, "~test direction~ utg%.6d , nothing after end\n", vu+1);
        // return 1;
        return 0;
    }else{
        fprintf(stderr, "~test direction~ utg%.6d , orphan\n", vu+1);
        return -1; 
    }
}

void hamt_utgutil_mark_graphdirection(ma_ug_t *ug, asg_t *sg){
    // purpose: traverse the unitig graph, mark direction of each tig
    //          (this function is meant to be used with ug DFS only, for the purpose of seeding in the correct direction)
    if (!ug->dir){  // ma_ug_gen uses calloc, so it's ok to test the pointer
        ug->dir = (int*)calloc(ug->u.n, sizeof(int));
    }
    memset(ug->dir, -1, sizeof(int)*ug->u.n);
    
    int base_dir = 0, cnt = 0;
    stacku32_t stack;
    stacku32_init(&stack);
    uint32_t vu, wu, v, w, v_;
    int32_t nv, nv_;
    asg_arc_t *av, *av_;

    stacku32_push(&stack, 0);
    ug->dir[0] = base_dir;
    uint8_t *color = (uint8_t*)calloc(ug->u.n, sizeof(uint8_t)); assert(color);

    uint32_t *v2vu = (uint32_t*)calloc(sg->n_seq*2, sizeof(uint32_t));
    uint8_t *sanmask = (uint8_t*)calloc(sg->n_seq*2, sizeof(uint8_t));
    for (uint32_t i=0; i<ug->u.n; i++){
        v2vu[ug->u.a[i].start] = i<<1 | 0;  // vertex to utg start/end status; similar to ma_ug_gen's mark
        v2vu[ug->u.a[i].end] = i<<1 | 1;
        sanmask[ug->u.a[i].start] = 1;  // just wanna make sure to prevent the case where i=0 and the vertex is a start...
        sanmask[ug->u.a[i].end] = 1;
    }

    while (1){  // main loop of traversal
        if (stacku32_pop(&stack, &vu)){
            fprintf(stderr, "[D::%s] pop, at utg%.6d\n", __func__, (int)vu+1);
            if (color[vu]==2){
                continue;
            }
            base_dir = ug->dir[vu];

            // examine the current node
            if (base_dir==0){
                v = ug->u.a[vu].end^1;
                v_ = ug->u.a[vu].start^1;
            }else{
                v = ug->u.a[vu].start^1;
                v_ = ug->u.a[vu].end^1;
            }
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            nv_ = asg_arc_n(sg, v_);
            av_ = asg_arc_a(sg, v_);
            if (nv==0 && nv_==0){  // is orphan
                color[vu] = 2;
                cnt++;
            }else{ // check child nodes and collect any white ones IN BOTH DIRECTION
                int is_exhausted = 1;

                // regular direction
                for (uint32_t i=0; i<nv; i++){  
                    w = av[i].v;
                    if (!sanmask[w]){continue;}
                    wu = v2vu[w];
                    if (color[wu>>1]!=0){continue;}
                    if (is_exhausted==1){
                        is_exhausted = 0;
                        stacku32_push(&stack, vu);  // put it back
                    }
                    if ((wu&1)==0){
                        ug->dir[wu>>1] = 0;
                    }else{
                        ug->dir[wu>>1] = 1;
                    }
                    color[wu>>1] = 1;
                    stacku32_push(&stack, wu>>1);
                    fprintf(stderr, "[D::%s]   push utg%.6d, dir is %d\n", __func__, (int)(wu>>1)+1, ug->dir[wu>>1]);
                }
                // reverse direction
                for (uint32_t i=0; i<nv_; i++){  
                    w = av_[i].v;
                    if (!sanmask[w]){continue;}
                    wu = v2vu[w];
                    if (color[wu>>1]!=0){continue;}
                    if (is_exhausted==1){
                        is_exhausted = 0;
                        stacku32_push(&stack, vu);  // put it back
                    }
                    if ((wu&1)==1){
                        ug->dir[wu>>1] = 0;
                    }else{
                        ug->dir[wu>>1] = 1;
                    }
                    color[wu>>1] = 1;
                    stacku32_push(&stack, wu>>1);
                    fprintf(stderr, "[D::%s]   push utg%.6d, dir is %d\n", __func__, (int)(wu>>1)+1, ug->dir[wu>>1]);
                }

                if (is_exhausted){
                    color[vu] = 2;
                    cnt++;
                }
            }
        }else{  // every time we get a new seed, it should be in a new subgraph
            if (cnt==ug->u.n){
                break;
            }
            base_dir = 0;
            for (vu=0; vu<ug->u.n; vu++){
                if (color[vu]==0){
                    break;
                }
                if (color[vu]==1){
                    fprintf(stderr, "a GRAY node.\n");
                    exit(1);
                }
            }
            if (vu==ug->u.n){
                fprintf(stderr, "abnormal traversal: can't find a new seed.\n");
                exit(1);
            }
            color[vu] = 1;
            ug->dir[vu] = base_dir;
            stacku32_push(&stack, vu);
            fprintf(stderr, "[D::%s] seed utg%.6d\n", __func__, (int)vu+1);
            
        } 
    }
    if (cnt==ug->u.n){
        fprintf(stderr, "[D::%s] finished traversal normally.\n", __func__);
    }else{
        fprintf(stderr, "[D::%s] abnormal traversal.\n", __func__);
        exit(1);
    }
    for (int i=0; i<ug->u.n; i++){
        fprintf(stderr, "utgdir - utg%.6d\t%d\n", i+1, ug->dir[i]);
        if (ug->dir[i]==-1){
            fprintf(stderr, "ERROR, utg direction not specified after traversal finished normally!\n");
            exit(1);
        }
    }

}


void hamt_utgutil_DFS(ma_ug_t *ug0, asg_t *sg0, uint64_t *ts_curr, uint64_t *ts_order, int *SCC_marking, int is_use_reverse){
    // NOTE: is only intended to be used in SCC routines; generic usage not guaranteed to be correct

    // original or transposed?
    ma_ug_t *ug;
    asg_t *sg;
    if (is_use_reverse){
        // ug = hamt_utgutil_transpose(ug0, sg0);
        // sg = ug->g;
        ug = ug0;  // don't change unitig graph
        sg = hamt_asggraph_util_gen_transpose(sg0);  // reverse all arcs
    }else{
        ug = ug0;
        sg = sg0;
    }

    // aux
    fprintf(stderr, "[D::%s] (entered).\n", __func__); fflush(stderr);
    double startTime = Get_T();
    int *color = (int*)calloc(ug->u.n, sizeof(int));  // 0 for white, 1 for gray, 2 for black
    stacku32_t DFSstack;
    stacku32_init(&DFSstack);
    uint32_t nb_node_finished = 0;
    uint64_t DFStime = 0;
    uint32_t nv;
    asg_arc_t *av;
    int SCC_label = 1;

    uint32_t v, w;
    uint32_t vu, wu;  // corresponding unitig index
    int dir;

    // creat mapping: asgNodeID to utgID
    uint8_t *sanmask = (uint8_t*)calloc(sg->n_seq*2, sizeof(uint8_t));
    uint32_t *v2vu = (uint32_t*)calloc(sg->n_seq*2, sizeof(uint32_t));
    uint8_t *utgdir = (uint8_t*)calloc(ug->u.n, sizeof(uint8_t));

    if (!is_use_reverse){
        for (uint32_t i=0; i<ug->u.n; i++){
            v2vu[ug->u.a[i].start] = i<<1 | 0;  // vertex to utg start/end status; similar to ma_ug_gen's mark
            v2vu[ug->u.a[i].end] = i<<1 | 1;
            sanmask[ug->u.a[i].start] = 1;
            sanmask[ug->u.a[i].end] = 1;
        }
    }else{
        // set utg direction
        for (int i=0; i<ug->u.n; i++){
            utgdir[((uint32_t)ts_order[i])>>1] = ((uint32_t)ts_order[i]) & 1;
        }
        for (uint32_t i=0; i<ug->u.n; i++){
            // because after asg transpose, end^1=>start will be start=>end^1
            //                              via av[i].v we will get utg.end^1 instead of utg.end
            v2vu[ug->u.a[i].start^1] = i<<1 | 1;  // if the target is the start of a unitig, that unitig is backward  
            v2vu[ug->u.a[i].end^1] = i<<1 | 0;  // if the target is the start of a unitig, that unitig is forward  
            sanmask[ug->u.a[i].start^1] = 1;
            sanmask[ug->u.a[i].end^1] = 1;
        }
    }

    // init
    if (ts_order){
        vu = (uint32_t)ts_order[ug->u.n-1];  // vu with direction
    }else{
        // mark all orphans
        uint32_t tmp=0;
        for (vu=0; vu<ug->u.n; vu++){  // here vu doesn't have direction
            dir = test_utg_direction_require_target(ug, sg, vu, sanmask);
            if (dir==-1){
                DFStime++;
                ts_curr[vu] = DFStime<<32 | vu<<1;
                if (dir==-1)
                    {fprintf(stderr, "ftA\tutg%.6d\t%d\tmarked-as-orphan\n", (int)vu+1, (int)DFStime);}  // finish time
                else
                    {fprintf(stderr, "ftA\tutg%.6d\t%d\tmarked-as-end\n", (int)vu+1, (int)DFStime);}  // finish time
                color[vu] = 2;
                nb_node_finished++;
            }else{
                tmp = vu<<1 | ug->dir[vu];
            }
        }
        vu = tmp; // vu with direction
        // vu = (uint32_t)0;
        // vu = (uint32_t) 4884<<1 | 0; // DEBUG
        // vu = (uint32_t)376<<1 | 1;  // DEBUG
        // vu = (uint32_t)949<<1 | 0;  // DEBUG
    }
    
    color[vu>>1] = 1;
    stacku32_push(&DFSstack, vu);
    if (SCC_marking){
        SCC_marking[vu>>1] = SCC_label;
    }

    
    // DFS main loop
    while (1){
        DFStime++;
        if (stacku32_pop(&DFSstack, &vu)){  // stack not empty, process the top one
            // vu has direction!

            fprintf(stderr, "pop, at utg%.6d, stack size is %d.\n", (int)(vu>>1)+1, DFSstack.n);
            if (color[vu>>1]==2){
                // this node has been colored black via other traversal
                DFStime--;
                continue;
            }

            // examine the current node
            if (!is_use_reverse){
                if (vu & 1){  // unitig of the reverse direction
                    v = ug->u.a[vu>>1].start^1; //////////UGH!!
                }else{
                    v = ug->u.a[vu>>1].end^1; ///////////UGH
                }
            }else{
                if (utgdir[vu>>1]){
                    v = ug->u.a[vu>>1].end;
                }else{
                    v = ug->u.a[vu>>1].start;  // unitig is forward, we search transpose-backward, so using the start vertex 
                }
            }
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            if (nv==0){  // is a leaf node, color it
                if (ts_curr) {
                    ts_curr[vu>>1] = DFStime<<32 | vu;
                    fprintf(stderr, "ftA\tutg%.6d\t%d\n", (int)(vu>>1)+1, (int)DFStime);  // finish time
                }
                color[vu>>1] = 2;
                nb_node_finished++;
                DFStime++;  //!!!
            }else{  // not a leaf node,
                int is_exhausted = 1;
                for (uint32_t i=0; i<nv; i++){  // check if there's any white child node
                    w = av[i].v;
                    if (!sanmask[w]){continue;}  // has arc, but the target isn't a valid utg start/end
                    if (color[v2vu[w]>>1]==0){
                        is_exhausted = 0;
                        break;
                    }
                }
                if (!is_exhausted){  // no, this node isn't finished
                    stacku32_push(&DFSstack, vu);  // put the unitig back
                    // push white child nodes
                    for (uint32_t i=0; i<nv; i++){
                        w = av[i].v;
                        fprintf(stderr, "   ...looking at node %.*s, utg%.6d\n", (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1), (int)(v2vu[w]>>1)+1);
                        if (!sanmask[w]){fprintf(stderr, "      (sanmask failed)"); continue;}  // has arc, but the target isn't a valid utg start/end
                        wu = v2vu[w];  // wu has direction
                        if ((ug->u.a[wu>>1].start>>1)==(ug->u.a[wu>>1].end>>1)){  // utg only has 1 read
                            if (!is_use_reverse){
                                if (ug->u.a[wu>>1].start==w){
                                    wu = wu>>1<<1 | 0;
                                }else{
                                    wu = wu>>1<<1 | 1;
                                }
                            }else{
                                if ((ug->u.a[wu>>1].start^1)==w){
                                    wu = wu>>1<<1 | 1;
                                }else{
                                    wu = wu>>1<<1 | 0;
                                }
                            }
                            // wu = wu>>1<<1 | (vu&1);  // use vu's direction
                        }
                        if (color[wu>>1]==0){
                            color[wu>>1] = 1;
                            stacku32_push(&DFSstack, wu);
                            DFStime++; // !!!!!!
                            if(SCC_marking) {SCC_marking[wu>>1] = SCC_label;}
                            fprintf(stderr, "  push, utg%.6d, direction is %d, target node %d\n", (int)((wu>>1)+1), wu&1, w>>1);
                            // debug
                            uint32_t xxx;
                            if (!is_use_reverse){
                                if (wu&1){
                                    xxx = ug->u.a[wu>>1].start^1;
                                }else{
                                    xxx = ug->u.a[wu>>1].end^1;
                                }
                            }else{
                                if (wu&1){
                                    xxx = ug->u.a[wu>>1].end;
                                }else{
                                    xxx = ug->u.a[wu>>1].start;
                                }
                            }
                            asg_arc_t *tmpav = asg_arc_a(sg, xxx);
                            for (int tmpi=0; tmpi<asg_arc_n(sg, xxx); tmpi++){
                                if (!sanmask[tmpav[tmpi].v]){continue;}
                                fprintf(stderr, "    (target: %d, direction %d)\n", (int)(v2vu[tmpav[tmpi].v]>>1)+1, (int)(v2vu[tmpav[tmpi].v]&1));
                            }
                        }
                    }
                }else{  // ok, this node is done (subtree all finished)
                    if (ts_curr){
                        ts_curr[vu>>1] = DFStime<<32 | vu;
                        fprintf(stderr, "ftA\tutg%.6d\t%d\n", (int)(vu>>1)+1, (int)DFStime);  // finish time
                    }
                    color[vu>>1] = 2;
                    nb_node_finished++;
                }
            }

        }else{  // stack is empty, get a new one or terminate
            if (nb_node_finished==ug->u.n){  // all done
                fprintf(stderr, "DFS all done, leaving\n");
                break;
            }else{  // get a new seed
                fprintf(stderr, "GET A NEW SEED\n");
                if (!ts_order){  // randomly get one seed
                    fprintf(stderr, "(random seed)\n");
                    uint32_t i;
                    for (i=0; i<ug->u.n; i++){
                        if (color[i]==0){
                            // vu = i<<1 | test_utg_direction_require_target(ug, sg, i, sanmask);
                            // vu = i<<1 | 0;
                            vu = i<<1 | ug->dir[i];
                            break;
                        }
                    }
                    if (i==ug->u.n){
                        // shouldn't happen
                        fprintf(stderr, "DEBUG WARNING, DFS random seed: leaving DFS with nb_node_finished=%d, total utg count=%d.\n", nb_node_finished, (int)ug->u.n);
                        break;
                    }
                }else{  // use the given order
                    fprintf(stderr, "(ordered seed)\n");
                    int i;
                    uint32_t j;
                    for (i=ug->u.n-1; i>=0; i--){
                        j = (uint32_t) ts_order[i];
                        if (color[j>>1]==0){
                            vu = j;
                            fprintf(stderr, "seed ftB\tutg%.6d\t%d\t%d\n", (int)(vu>>1)+1, (int)vu&1,(int)(ts_order[i]>>32));  // finish time
                            break;
                        }
                    }
                    if (i==-1){
                        // shouldn't happen
                        fprintf(stderr, "DEBUG WARNING, DFS ordered seed: leaving DFS with nb_node_finished=%d, total utg count=%d.\n", nb_node_finished, (int)ug->u.n);
                        break;
                    }
                }
                // update status
                fprintf(stderr, "seeding, at %d, stack size is %d.\n", (int)vu>>1, DFSstack.n);
                color[vu>>1] = 1;
                stacku32_push(&DFSstack, vu);
                SCC_label++;
                if (SCC_marking){
                    SCC_marking[vu>>1] = SCC_label;
                    fprintf(stderr, "CURRENT SCC LABEL is: %d", SCC_label);
                }

            }

        }

    }

    // clean up
    if (is_use_reverse){
        asg_destroy(sg);
        // ma_ug_destroy(ug);  // free the transposed graph
    }
    free(v2vu);
    free(sanmask);
    stacku32_destroy(&DFSstack);
    free(color);

    fprintf(stderr, "[D::%s] (exited).\n", __func__); fflush(stderr);
    fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);

}

int *hamt_utg_SCC(ma_ug_t *ug, asg_t *sg){
    // aux
    fprintf(stderr, "[D::%s] (entered).\n", __func__); fflush(stderr);
    double startTime = Get_T();
    asg_cleanup(sg);

    uint64_t *ts1 = (uint64_t*)calloc(ug->u.n, sizeof(uint64_t));  // upper 32 bits: finish timestamp, lower 32 bits: node ID
    int *SCC_labels = (int*)calloc(ug->u.n, sizeof(int));

    // 1st DFS
    hamt_utgutil_DFS(ug, sg, ts1, NULL, NULL, 0);

    // get the ordering (finish timestamps as of in 1st DFS)
    radix_sort_ovhamt64(ts1, ts1+ug->u.n);
    
    // 2nd DFS
    hamt_utgutil_DFS(ug, sg, NULL, ts1, SCC_labels, 1);

    free(ts1);
    if (VERBOSE>=1){
        fprintf(stderr, "[D::%s] (exited).\n", __func__);
        fprintf(stderr, "[M::%s] took %0.2fs\n\n", __func__, Get_T()-startTime);
    }
    return SCC_labels;
}


void hamt_utg_SCCcovCut(asg_t *sg,
                        const ma_sub_t* coverage_cut,   // for utg coverage
                        ma_hit_t_alloc* sources, R_to_U* ruIndex)  // for utg coverage
{
    // (do cutting on the asg, then regenerate ug)
    // mark strongly connected components at unitig graph level, 
    // and cut the arcs with differential coverages that connect two SCCs

    asg_cleanup(sg);
    ma_ug_t *ug = ma_ug_gen(sg);
    hamt_utgutil_mark_graphdirection(ug, sg);  // mark unitig directions, for DFS to pick the correct seed direction

    fprintf(stderr, "[debug::%s] n_utg is %d.", __func__, (int)ug->u.n);

    // aux
    int n_reduced = 0;
    uint32_t v, w, ww, n_vtx = sg->n_seq*2, nv, nv_;
    uint32_t vu, wu;
    asg_arc_t *av, *av_;
    int SCC_label_v, SCC_label_w;
    int utgcov_v, utgcov_w;
    uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));  // for utg coverage
    int dummy;  // for breaking-out-of-if-clause


    // SCC labeling
    int *SCC_labels = hamt_utg_SCC(ug, sg);
    
    if (VERBOSE>=1){
        // debug
        int max_SCC_label = 1;
        for (int i=0; i<ug->u.n; i++){
            if (SCC_labels[i]>max_SCC_label){
                max_SCC_label = SCC_labels[i];
            }
        }
        fprintf(stderr, "[D::%s] asg has %d strongly conneccted components.\n", __func__, max_SCC_label);
        write_debug_SCC_labeling_ug(ug, SCC_labels, asm_opt.output_file_name);

        char tmp[100];
        sprintf(tmp, "%s.SCCug.gfa", asm_opt.output_file_name);
        FILE *fp = fopen(tmp, "w");
        ma_ug_print2_lite(ug, &R_INF, sg, 0, "utg", fp);
        fclose(fp);
    }
    int *SCC_labels_vertices = (int*)malloc(sizeof(int) * sg->n_seq);
    memset(SCC_labels_vertices, -1, sizeof(int) * sg->n_seq);  // -1 for sanity check
    uint32_t *v2vu = (uint32_t*)malloc(sizeof(uint32_t) * sg->n_seq);
    for (vu=0; vu<ug->u.n; vu++){
        SCC_labels_vertices[ug->u.a[vu].start>>1] = SCC_labels[vu];
        SCC_labels_vertices[ug->u.a[vu].end>>1] = SCC_labels[vu];
        v2vu[ug->u.a[vu].start>>1] = vu;  // v2vu here is just for helping to check utg coverage, we don't need direction info
        v2vu[ug->u.a[vu].end>>1] = vu;
    }

    // SCC-aware coverge-based cutting
    uint32_t placeholder[2];
    for (vu=0; vu<ug->u.n; vu++){
        fprintf(stderr, "COVCUT: at utg #%d\n", vu+1);
        SCC_label_v = SCC_labels[vu];
        utgcov_v = (int) get_ug_coverage(&ug->u.a[vu], sg, coverage_cut, sources, ruIndex, primary_flag);
        placeholder[0] = ug->u.a[vu].start^1;  // test the predecessor
        placeholder[1] = ug->u.a[vu].end^1;  // test the target

        int idx_vu = 0;
        while (idx_vu<2){
            dummy = 1;

            while(dummy){
                v = placeholder[idx_vu];

                nv = asg_arc_n(sg, v);
                av = asg_arc_a(sg, v);

                if (nv==0){
                    fprintf(stderr, "> utg%.6d no target\n", vu+1);
                    fprintf(stderr, "  (bcs nv)\n");
                    break;
                }
                for (int idx_target=0; idx_target<nv; idx_target++){  // check each target
                    w = av[idx_target].v;
                    fprintf(stderr, "> utg%.6d ==> utg%.6d\n", vu+1, v2vu[w>>1]+1);

                    SCC_label_w = SCC_labels_vertices[w>>1];
                    if (SCC_label_v==SCC_label_w){
                        fprintf(stderr, "  (bcs same SCC label)\n");
                        // break;  // dummy loop
                        continue;
                    }

                    assert(SCC_label_w!=-1);
                    
                    // require that w has a predecessor that shares the same SCClabel with w (i.e. testing the loop)
                    nv_ = asg_arc_n(sg, w^1);
                    // if (nv_!=2){   // TODO: allow more
                    if (nv_==0){
                        fprintf(stderr, "  (bcs nv_)\n");
                        // break;  // dummy loop
                        continue;
                    }
                    av_ = asg_arc_a(sg, w^1);
                    int looks_like_a_loop = 0;
                    for (int idx_rev_target=0; idx_rev_target<nv_; idx_rev_target++){
                        ww = av_[0].v;
                        if ((ww>>1)==(v>>1)){ 
                            ww = av_[1].v;
                            looks_like_a_loop = 1;
                            break;
                        }
                    }
                    if (!looks_like_a_loop){
                        fprintf(stderr, "  (bcs secondary SCC label)\n");
                        // break; // dummy loop
                        continue;
                    }

                    // check coverage diff
                    wu = v2vu[w>>1];
                    utgcov_w = (int) get_ug_coverage(&ug->u.a[wu], sg, coverage_cut, sources, ruIndex, primary_flag);
                    int please_cut = 0;
                    if (utgcov_w>=10){
                        if (utgcov_w<20){
                            if ((float)utgcov_w/utgcov_v>3 || (utgcov_w-utgcov_v>10))
                                please_cut = 1;
                        }else{
                            if ((float)utgcov_w/utgcov_v>3)
                                please_cut = 1;
                        }
                    }
                    if (!please_cut){  // don't cut
                        fprintf(stderr, "  (bcs target coverage)\n");
                        // break; // dummy loop
                        continue;
                    }
            
                    // cut the arc
                    asg_arc_del(sg, v, w, 1);
                    asg_arc_del(sg, w^1, v^1, 1);
                    n_reduced++;
                    fprintf(stderr, "  CUT!\n");
                    if (VERBOSE>=1){
                        fprintf(stderr, "[Vdebug::%s] SCC-aware cov-based cut (%d vs %d).\n", __func__, utgcov_v, utgcov_w);
                    }
                }

                dummy=0;  // for leaving the dummy while loop
            }
            idx_vu++;  // go check the other end
        }
    }

    if (n_reduced){
        fprintf(stderr, "[M::%s] dropped %d links.\n", __func__, n_reduced);
        asg_cleanup(sg);
    }


    // // redo ug
    // // (TODO: is there a better way to update unitig graph? (note: the current code alters the working sg, not the ug->g))
    // hamt_utg_redoug(ug, sg);
    // ma_ug_destroy(ug);
    free(primary_flag);
    free(SCC_labels);
    free(SCC_labels_vertices);
    free(v2vu);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////// experimental SCC, last try before i give up ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void hamt_mark_usg_v2vu_nodirection(ma_ug_t *ug, asg_t *usg){
    if (usg->v2vu){
        free(usg->v2vu);
    }
    usg->v2vu = (uint32_t*)calloc(usg->n_seq, sizeof(uint32_t));  assert(usg->v2vu);  // if you segfault here: did you forgot to >>1
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        usg->v2vu[ug->u.a[vu].start>>1] = vu;
        usg->v2vu[ug->u.a[vu].end>>1] = vu;
    }
}

// int hamt_check_usg_weirdarc(asg_t *usg, uint32_t v, uint32_t w){
//     assert(usg->v2vu);
//     uint32_t u0, u1, u2;
//     u0 = usg->v2vu[v>>1];
//     u1 = usg->v2vu[w>>1];
//     if (u0==u1){
//         return 0;
//     }
//     asg_arc_t *av = asg_arc_a(usg, w^1);
//     uint32_t nv = asg_arc_n(usg, w^1);
//     if (nv==0){
//         return 0;
//     }
//     for (int i=0; i<nv; i++){
//         u2 = usg->v2vu[av[i].v>>1];
//         if (u1!=u2 && u2!=u0){
//             return 1;
//         }
//     }
//     return 0;
// }

void hamt_mark_usg_direction(asg_t *usg){
    // chose arbitrary traversal direction for each node, aka block the complementary arcs
    assert(usg->v2vu);
    if (usg->whitelist){  // note: asg_init uses calloc, this is safe
        free(usg->whitelist);
        free(usg->blacklist);
    }
    usg->whitelist = (uint8_t*)calloc(usg->n_seq*2, sizeof(uint8_t)); assert(usg->whitelist);
    usg->blacklist = (uint8_t*)calloc(usg->n_seq*2, sizeof(uint8_t)); assert(usg->blacklist);
    
    // aux
    uint8_t *color = (uint8_t*)calloc(usg->n_seq*2, sizeof(uint8_t)); assert(color); 
    uint8_t *color_node = (uint8_t*)calloc(usg->n_seq, sizeof(uint8_t)); assert(color); 
    stacku32_t stack;
    stacku32_init(&stack);

    // first seed
    uint32_t v;
    for (v=0; v<usg->n_seq*2; v++){
        if (asg_arc_n(usg, v)>0){
            break;
        }
    }
    stacku32_push(&stack, v);
    color[v] = 1;
    color_node[v>>1]++;
    usg->whitelist[v] = 1;

    uint32_t nv, w;
    asg_arc_t *av;
    while (1){  // main loop of DFS
        if (stacku32_pop(&stack, &v)){
            if (color[v]==2){
                continue;
            }
            uint32_t i;
            int is_exhausted = 1;
            int is_weird;
            // search forward direction
            nv = asg_arc_n(usg, v);
            av = asg_arc_a(usg, v);
            for (i=0; i<nv; i++){
                w = av[i].v;
                if (color[w]==0){
                    // ////// check if this is a weird direction
                    // is_weird = hamt_check_usg_weirdarc(usg, v, w);
                    // if (is_weird){
                    //     usg->blacklist[w] = 1;
                    //     fprintf(stderr, "[D::%s] ignoring an arc between utg%.6d and utg%.6d.\n", __func__, (int)usg->v2vu[v>>1]+1, (int)usg->v2vu[w>>1]+1);
                    //     fprintf(stderr, "        (read %.*s to read %.*s.)\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1));
                    //     continue;  // pretend that we didn't see this node
                    // }
                    if (color[w^1]!=0){
                        usg->blacklist[w] = 1;
                    }
                    // ///// end of check
                    if (is_exhausted){
                        is_exhausted = 0;
                        stacku32_push(&stack, v);  // put back the current node
                    }
                    stacku32_push(&stack, w);
                    color[w] = 1;
                    color_node[w>>1]++;
                    usg->whitelist[w] = 1;
                }
            }
            // search backward direction (note vertices are added in the forward direction)
            nv = asg_arc_n(usg, v^1);
            av = asg_arc_a(usg, v^1);
            for (i=0; i<nv; i++){
                w = av[i].v^1;
                if (color[w]==0){
                    // ////// check if this is a weird direction
                    // is_weird = hamt_check_usg_weirdarc(usg, v, w);
                    // if (is_weird){
                    //     usg->blacklist[w] = 1;
                    //     fprintf(stderr, "[D::%s] ignoring an arc between utg%.6d and utg%.6d.\n", __func__, (int)usg->v2vu[v>>1]+1, (int)usg->v2vu[w>>1]+1);
                    //     fprintf(stderr, "        (read %.*s to read %.*s.)\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1));
                    //     continue;  // pretend that we didn't see this node
                    // }
                    // ///// end of check
                    if (is_exhausted){
                        is_exhausted = 0;
                        stacku32_push(&stack, v);  // put back the current node
                    }
                    stacku32_push(&stack, w);
                    color[w] = 1;
                    color_node[w>>1]++;
                    usg->whitelist[w] = 1;
                }
            }
        }else{  // finished one subgraph
            uint32_t i=0;
            int is_finished = 1;
            for (i=0; i<usg->n_seq; i++){
                if (color_node[i]==0){
                    if (asg_arc_n(usg, i<<1)>0){
                        v = i<<1;
                        is_finished = 0;
                        break;
                    }else if (asg_arc_n(usg, i<<1 | 1)>0){
                        v = i<<1 | 1;
                        is_finished = 0;
                        break;
                    }
                }
            }
            if (i==usg->n_seq){
                if (!is_finished){
                    fprintf(stderr, "[E::%s] abnormal finish. \n", __func__);
                    exit(1);
                }
                fprintf(stderr, "[D::%s] normal finish. \n", __func__);
                break;  // finish
            }
            stacku32_push(&stack, v);
            color[v] = 1;
            color_node[v>>1]++;
            usg->whitelist[v] = 1;
        }
    }
    stacku32_destroy(&stack);
    free(color);
    free(color_node);
}

asg_t *hamt_abstract_ug_to_sg(asg_t *sg0, ma_ug_t *ug){
    // (meant to be used with ug SCC.)
    // generate, sort and index a string graph that's only contains start/end nodes of the unitig graph,
    // and their corresponding linkages (including dummy bidirection ones between start and end of a utg)
    asg_t *g = asg_init();
    uint32_t i, nv;
    asg_arc_t *ap, *aq;
    uint32_t *v2vu = (uint32_t*)calloc(sg0->n_seq, sizeof(uint32_t));
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        v2vu[ug->u.a[vu].start>>1] = vu;
        v2vu[ug->u.a[vu].end>>1] = vu;
    }

    uint8_t *mark = (uint8_t*)calloc(sg0->n_seq*2, sizeof(uint8_t));
    for (i=0; i<ug->u.n; i++){
        mark[ug->u.a[i].start] = 1;
        mark[ug->u.a[i].end] = 1;
        mark[ug->u.a[i].start^1] = 1;
        mark[ug->u.a[i].end^1] = 1;
        if ((ug->u.a[i].start>>1)!=(ug->u.a[i].end>>1)){
            g->usg_n_seq+=2;
        }else{
            g->usg_n_seq+=1;
        }
    }
    g->n_seq = sg0->n_seq;  // note: it's sparse but we do need the whole bucket to keep IDs consistent (with min effort)

    // add arcs between unitigs
    for (i=0; i<sg0->n_arc; i++){
        ap = &sg0->arc[i];
        if (mark[ap->ul>>32]>0 && mark[ap->v]>0){
            // if this happens to be a internal arc for a 2-read unitig, skip it
            if (v2vu[ap->ul>>33]==v2vu[ap->v>>1]){
                continue;
            }

            // (the other direction will be added during other iters)
            aq = asg_arc_pushp(g);
            aq->ol = ap->ol;
            aq->del = 0;
            aq->ul = ap->ul;
            aq->v = ap->v;
        }
    }

    uint32_t base_i = i;
    uint32_t start, end;
    // add dummy arcs between start and end of a unitig
    for (i=0; i<ug->u.n; i++){
        start = ug->u.a[i].start;
        end = ug->u.a[i].end;
        // // skip 1-read unitig
        // if ((start>>1)==(end>>1)){
        //     continue;
        // }

        // forward, direction 1
        aq = asg_arc_pushp(g);
        aq->ol = 10; // arbitrary
        aq->del = 0;
        aq->ul = ((uint64_t)start)<<32 | 1;
        aq->v = end^1;

        // forward, direction 2
        aq = asg_arc_pushp(g);
        aq->ol = 10; // arbitrary
        aq->del = 0;
        aq->ul = ((uint64_t)end)<<32 | 1;
        aq->v = start^1;

        // reverse, direction 1
        aq = asg_arc_pushp(g);
        aq->ol = 10; // arbitrary
        aq->del = 0;
        aq->ul = ((uint64_t)(end^1))<<32 | 11;
        aq->v = start;

        // reverse, direction 2
        aq = asg_arc_pushp(g);
        aq->ol = 10; // arbitrary
        aq->del = 0;
        aq->ul = ((uint64_t)(start^1))<<32 | 1;
        aq->v = end;
    }
    
    // TODO: check symmetry

    asg_arc_sort(g);
    g->is_srt = 1;
    fprintf(stderr, "usg n_seq is %d (actual node number: %d), n_arc is %d.\n", (int)g->n_seq, (int)g->usg_n_seq, (int)g->n_arc);
    asg_arc_index(g);
    
    free(v2vu);

    return g;
}


void hamt_usg_DFS(ma_ug_t *ug, asg_t *usg0, uint64_t *ts_curr, uint64_t *ts_order, int *SCC_marking, int is_use_reverse){
    // transpose?
    asg_t *usg;
    if (is_use_reverse){
        usg = hamt_asggraph_util_gen_transpose(usg0);  // reverse all arcs
        usg->whitelist = usg0->whitelist;
    }else{
        usg = usg0;
    }

    // aux
    stacku32_t stack;
    stacku32_init(&stack);
    uint8_t *color = (uint8_t*)calloc(usg->n_seq*2, sizeof(uint8_t));
    int SCC_label = 1;
    uint64_t DFStime = 0;
    int nb_node_finished = 0;

    uint32_t *v2vu = (uint32_t*)calloc(usg->n_seq, sizeof(uint32_t));
    for (uint32_t vu=0; vu<ug->u.n; vu++){
        v2vu[ug->u.a[vu].start>>1] = vu;
        v2vu[ug->u.a[vu].end>>1] = vu;
    }

    
    // handles
    uint32_t nv, v, w;
    asg_arc_t *av;

    // init finishing time
    if (ts_curr){
        for (uint32_t i=0; i<usg->n_seq*2; i++){
            ts_curr[i] = (uint64_t)i;
        }
    }

    // init DFS
    if (!is_use_reverse){  // pick a random seed
        assert(ts_curr);
        for (v=0; v<usg->n_seq*2; v++){
            if (!usg->whitelist[v]){
                continue;
            }
            if (asg_arc_n(usg, v)>0){
                break;
            }
        }
        if (v==usg->n_seq*2){
            fprintf(stderr, "[E::%s] can't get the 1st random seed.\n", __func__);
            exit(1);
        }
    }else{  // pick a ordered seed
        assert(ts_order);
        v = (uint32_t)ts_order[usg->n_seq*2-1];
        if (SCC_marking){
            SCC_marking[v] = SCC_label;
        }
    }
    stacku32_push(&stack, v);
    DFStime++;
    color[v] = 1;

    // DFS main
    while (1){
        DFStime++;
        if (stacku32_pop(&stack, &v)){
            assert(usg->whitelist[v]);
            fprintf(stderr, "pop utg%.6d ( %.*s dir %d)\n", (int)v2vu[v>>1]+1, (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)v&1);
            if (color[v]==2){
                DFStime--;
                continue;
            }
            // examine the current node
            nv = asg_arc_n(usg, v);
            av = asg_arc_a(usg, v);
            if (nv==0){  // leaf node, color it
                if(ts_curr){
                    ts_curr[v] = DFStime<<32 | v;
                }
                color[v] = 2;
                nb_node_finished++;
            }else{
                int is_exhausted = 1;
                for (uint32_t i=0; i<nv; i++){
                    w = av[i].v;
                    if (!usg->whitelist[w]){continue;}
                    if (color[w]==0){  // put the current node back
                        if (is_exhausted){
                            is_exhausted = 0;
                            stacku32_push(&stack, v);
                            fprintf(stderr, "    put back utg%.6d ( %.*s dir %d)\n", (int)v2vu[v>>1]+1, (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)v&1);
                        }
                        color[w] = 1;
                        DFStime++;
                        if (SCC_marking){
                            SCC_marking[w] = SCC_label;
                            fprintf(stderr, "SCCdebug\tutg%.6d\t%.*s\t%d\n", (int)v2vu[w>>1]+1, (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1), SCC_label);
                        }
                        stacku32_push(&stack, w);
                        fprintf(stderr, "    push utg%.6d ( %.*s dir %d)\n", (int)v2vu[w>>1]+1, (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1), (int)w&1);
                    }
                }
                if (is_exhausted){
                    color[v] = 2;
                    nb_node_finished++;
                    if (ts_curr){
                        ts_curr[v] = DFStime<<32 | v;
                    }
                }
            }
        }else{
            fprintf(stderr, "SEED\n");
            if (!is_use_reverse){  // random seed
                uint32_t i;
                for (i=0; i<usg->n_seq*2; i++){
                    if (color[i]==0 && usg->whitelist[i]){
                        break;
                    }
                }
                if (i==usg->n_seq*2){
                    // fprintf(stderr, "[E::%s] random seeding failed while there's nodes left. (%d / %d)\n", __func__, (int)nb_node_finished, (int)(usg->n_seq*2));
                    // exit(1);
                    break;
                }
                color[i] = 1;
                stacku32_push(&stack, i);
                v = i;
                fprintf(stderr, "    seeded utg%.6d ( %.*s dir %d)\n", (int)v2vu[v>>1]+1, (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)v&1);
            }else{  // ordered seed
                SCC_label++;
                int i;
                for (i=usg->n_seq*2-1; i>=0; i--){
                    if ((ts_order[i]>>32)==0){
                        continue;
                    }
                    if (!usg->whitelist[(uint32_t)ts_order[i]]){
                        continue;
                    }
                    if (color[(uint32_t)ts_order[i]]==0){
                        v = (uint32_t)ts_order[i];
                        break;
                    }
                }
                if (i==-1){
                    // fprintf(stderr, "[E::%s] ordered seeding failed while there's nodes left. (%d / %d)\n", __func__, (int)nb_node_finished, (int)(usg->n_seq*2));
                    // exit(1);
                    break;
                }
                color[v] = 1;
                stacku32_push(&stack, v);
                SCC_marking[v] = SCC_label;
                fprintf(stderr, "SCCdebug\tutg%.6d\t%.*s\t%d\n", (int)v2vu[v>>1]+1, (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), SCC_label);
                fprintf(stderr, "    seeded utg%.6d ( %.*s dir %d)\n", (int)v2vu[v>>1]+1, (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1), (int)v&1);
                fprintf(stderr, "    current SCC label: %d.\n", SCC_label);
            }
        }
    }

    if (is_use_reverse){
        asg_destroy(usg);
    }
    free(color);
    stacku32_destroy(&stack);
    fprintf(stderr, "finished DFS.\n\n\n");
    free(v2vu);
}


int *hamt_usg_SCC(ma_ug_t *ug, asg_t *usg){
    // aux
    fprintf(stderr, "[D::%s] entered.\n", __func__);
    uint64_t *ts1 = (uint64_t*)calloc(usg->n_seq*2, sizeof(uint64_t));
    int *SCC_labels = (int*)calloc(usg->n_seq*2, sizeof(int));
    assert(usg->whitelist);
    int san = 0;
    for (uint32_t i=0; i<usg->n_seq*2; i++){
        if (usg->whitelist[i])
            san++;
    }
    fprintf(stderr, "items in the whitelist: %d.\n", san);

    hamt_usg_DFS(ug, usg, ts1, NULL, NULL, 0);  // 1st DFS
    radix_sort_ovhamt64(ts1, ts1+usg->n_seq*2);  // sort by finishing time (ID packed)
    hamt_usg_DFS(ug, usg, NULL, ts1, SCC_labels, 1);  // 2nd DFS, transposed

    free(ts1);
    return SCC_labels;
}


void hamt_usg_SCCcovCut(asg_t *sg0,
                        const ma_sub_t* coverage_cut,   // for utg coverage
                        ma_hit_t_alloc* sources, R_to_U* ruIndex)  // for utg coverage
{
    asg_cleanup(sg0);
    ma_ug_t *ug0 = ma_ug_gen(sg0);
    asg_t *usg = hamt_abstract_ug_to_sg(sg0, ug0);
    hamt_mark_usg_v2vu_nodirection(ug0, usg);
    hamt_mark_usg_direction(usg);  // note: need v2vu(nodir) marked
    

    // debug
    char tmp[100];
    sprintf(tmp, "%s.SCCug.gfa", asm_opt.output_file_name);
    FILE *fp = fopen(tmp, "w");
    ma_ug_print2_lite(ug0, &R_INF, sg0, 0, "utg", fp);
    fclose(fp);
    // end of debug

    int *SCC_labels = hamt_usg_SCC(ug0, usg);

    int *SCC_labels_node = (int*)calloc(usg->n_seq, sizeof(int));
    int *coverage_on_usg = (int*)calloc(usg->n_seq, sizeof(int));
    for (uint32_t i=0; i<usg->n_seq*2; i++){  // TODO: merge this into SCC
        if (!usg->whitelist[i]){continue;}
        if (SCC_labels_node[i>>1]==0){
            SCC_labels_node[i>>1] = SCC_labels[i];
        }else{
            if (SCC_labels[i]!=0)
                // assert(SCC_labels_node[i>>1]==SCC_labels[i]);
                fprintf(stderr, "SCCwarning: inconsistent label? read is %.*s, node label %d, vertex label %d\n", (int)Get_NAME_LENGTH(R_INF, i>>1), Get_NAME(R_INF, i>>1), SCC_labels_node[i>>1], SCC_labels[i]);
        }
    }
    write_debug_SCC_labeling2(ug0, usg, SCC_labels, asm_opt.output_file_name);

    // get coverage info
    uint32_t vu, v;
    int label_v, utgcov_v, utgcov_w;
    uint8_t* primary_flag = (uint8_t*)calloc(sg0->n_seq, sizeof(uint8_t));  // for utg coverage
    int dummy;
    for (vu=0; vu<ug0->u.n; vu++){
        utgcov_v = (int) get_ug_coverage(&ug0->u.a[vu], sg0, coverage_cut, sources, ruIndex, primary_flag);
        v = ug0->u.a[vu].start>>1;
        coverage_on_usg[v] = utgcov_v;
        v = ug0->u.a[vu].end>>1;
        coverage_on_usg[v] = utgcov_v;
    }

    // cut arcs on the string graph, using the usg
    uint32_t s, e;
    asg_arc_t *ap;
    int nb_cut = 0;
    int yescut = 0;
    for (uint32_t i=0; i<usg->n_arc; i++){
        s = usg->arc[i].ul>>32;
        e = usg->arc[i].v;

        // check label diff
        if (SCC_labels_node[s>>1]==SCC_labels_node[e>>1]){
            continue;
        }

        // check unitig coverage diff
        yescut = 0;
        utgcov_v = coverage_on_usg[s>>1];
        utgcov_w = coverage_on_usg[e>>1];
        if (utgcov_w>=10){
            if (utgcov_w<20){
                if ((float)utgcov_w/utgcov_v>3 || (utgcov_w-utgcov_v>10))
                    yescut = 1;
            }else{
                if ((float)utgcov_w/utgcov_v>3)
                    yescut = 1;
            }
        }

        // cut
        if (yescut){
            nb_cut++;
            asg_arc_del(sg0, s, e, 1);
            asg_arc_del(sg0, e^1, s^1, 1);
        }
        
    }

    free(usg->whitelist);
    free(usg->blacklist);
    asg_destroy(usg);
    free(SCC_labels);
    free(SCC_labels_node);
    // free(SCC_labels_utg);
    free(coverage_on_usg);
    free(primary_flag);

    // clean graph
    if (nb_cut){
        asg_cleanup(sg0);
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////              ok NO.                         ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////