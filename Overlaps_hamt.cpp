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

        for (i=0; i<nv; i++)
            if (!av[i].del) 
                {v = av[0].v; 
                acc += (uint32_t) av[0].ul;
                break;}
        if (hamt_asgarc_util_countPre(sg, v, 0, 0)!=1) return 0;  // backward bifur 
    }
    return 1;
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
        a->ul = ((uint64_t)a0->v)<<32 | (uint64_t)((uint32_t)a0->ul);
        a->v = a0->ul>>32;
    }
    
    asg_arc_sort(g);
    g->is_srt = 1;
    asg_arc_index(g);

    return g;

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

int hamt_ugarc_detect_circle(ma_ug_t *ug, asg_t *sg, uint32_t vu0, uint32_t wu0, uint32_t *v2vu, uint8_t *sanmask){ // vu0 has direction; v2vu has directions
    // DFS on ug, test if vu is the start of a circle
    // return 1 if yes, 0 otherwise

    uint32_t vu=wu0, nv, v, w, wu;
    asg_arc_t *av;
    uint8_t *color = (uint8_t*)calloc(ug->u.n*2, 1);
    stacku32_t stack;
    stacku32_init(&stack);
    stacku32_push(&stack, vu);
    color[vu] = 1;

    int ret = 0;

    // fprintf(stderr, ">checking utg%.6d\n", (int)(vu>>1)+1);
    while (1){  // DFS main loop
        if (stacku32_pop(&stack, &vu)){
            // fprintf(stderr, "    pop utg%.6d (dir %d)\n", (int)(vu>>1)+1, (int)(vu&1));
            if (color[vu]==2){
                continue;
            }

            if ((vu&1)==0){
                v = ug->u.a[vu>>1].end^1;
            }else{
                v = ug->u.a[vu>>1].start^1;
            }
            // fprintf(stderr, "        handle is %.*s\n", (int)Get_NAME_LENGTH(R_INF, v>>1), Get_NAME(R_INF, v>>1));
            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            if (nv==0){
                color[vu] = 2;
                continue;
            }
            int is_exhausted = 1;
            for (int i=0; i<(int)nv; i++){
                w = av[i].v;
                wu = v2vu[w];

                if (wu==wu0){  // found the circle
                    ret = 1;  // don't exit, we need to empty the stack to make sure the handle tig isn't in
                }
                if (wu==vu0){  // the circle shall not include the handle tig
                    // fprintf(stderr, "HANDLE IN CIRCLE.\n");
                    free(color);
                    stacku32_destroy(&stack);
                    return 0;
                }
                if (!sanmask[w]){  // not a start/end of a unitig
                    continue;
                }

                if ((ug->u.a[wu>>1].start>>1)==(ug->u.a[wu>>1].end>>1)){  // utg only has 1 read
                    if (ug->u.a[wu>>1].start==w){
                        wu = wu>>1<<1 | 0;
                    }else{
                        wu = wu>>1<<1 | 1;
                    }
                }


                if (color[wu]==0){
                    if (is_exhausted){
                        is_exhausted = 0;
                        stacku32_push(&stack, vu); // put the current node back
                        // fprintf(stderr, "    put back utg%.6d\n", (int)(vu>>1)+1);
                    }
                    color[wu] = 1;
                    stacku32_push(&stack, wu);
                    // fprintf(stderr, "    pushed utg%.6d (dir %d)\n", (int)(wu>>1)+1, (int)(wu&1));
                    // fprintf(stderr, "        handle is %.*s\n", (int)Get_NAME_LENGTH(R_INF, w>>1), Get_NAME(R_INF, w>>1));
                }
            }
            if (is_exhausted){
                color[vu] = 2;
            }
        }else{
            break;
        }
    }
    free(color);
    stacku32_destroy(&stack);
    return ret;
}


void hamt_ugarc_covcut_danglingCircle(asg_t *sg,
                                      const ma_sub_t* coverage_cut,   // for utg coverage
                                      ma_hit_t_alloc* sources, R_to_U* ruIndex)  // for utg coverage
{
    ma_ug_t *ug = ma_ug_gen(sg);
    uint32_t vu, vu0, wu, nv, v, w;
    uint32_t cov1, cov2, covtmp;
    asg_arc_t *av;
    uint32_t *v2vu_directed = hamt_ugarc_util_v2vu_yesdirection(ug, sg);
    uint8_t *sanmask = (uint8_t*)calloc(sg->n_seq*2, sizeof(uint8_t));  // just in case a start ID is zero
    for (uint32_t i=0; i<ug->u.n; i++){
        sanmask[ug->u.a[i].start] = 1;
        sanmask[ug->u.a[i].end] = 1;
    }
    uint8_t* primary_flag = (uint8_t*)calloc(sg->n_seq, sizeof(uint8_t));  // for utg coverage

    if (VERBOSE>=1){
        char tmp[100];
        sprintf(tmp, "%s.SCCug.gfa", asm_opt.output_file_name);
        FILE *fp = fopen(tmp, "w");
        ma_ug_print2_lite(ug, &R_INF, sg, 0, "utg", fp);
        fclose(fp);
    }

    int ret;
    int nb_cut =0;
    for (vu=0; vu<ug->u.n; vu++){
        cov1 = get_ug_coverage(&ug->u.a[vu], sg, coverage_cut, sources, ruIndex, primary_flag);

        for (uint32_t dir=0; dir<=1; dir++){  // try both directions
            if (VERBOSE>=1) {fprintf(stderr, "@ utg%.6d, dir %d\n", vu+1, (int)(dir));}
            if (dir==0){
                v = ug->u.a[vu].end^1;
            }else{
                v = ug->u.a[vu].start^1;
            }

            nv = asg_arc_n(sg, v);
            av = asg_arc_a(sg, v);
            if (nv<1){
                if (VERBOSE>=1) {fprintf(stderr, "  (nv=0)\n");}
                continue;
            }
            for (int i=0; i<(int)nv; i++){
                w = av[i].v;
                if (!sanmask[w]){
                    if (VERBOSE>=1) {fprintf(stderr, "  (w is not unitig start/end.)\n");}
                    continue;  // not a unitig start/end
                }
                if (asg_arc_n(sg, w^1)<2){  // early termination, the target can't be a circle
                    if (VERBOSE>=1) {fprintf(stderr, "  (w has less than 2 backward arc.)\n");}
                    continue;
                }
                wu = v2vu_directed[w];
                if (VERBOSE>=1) {fprintf(stderr, "   [target is utg%.6d ]\n", (int)(wu>>1)+1);}
                cov2 = get_ug_coverage(&ug->u.a[wu>>1], sg, coverage_cut, sources, ruIndex, primary_flag);

                if (cov1>cov2){
                    covtmp = cov1;
                    cov1 = cov2;
                    cov2 = covtmp;
                }

                if (cov2<10){  // early termination, target tig has low coverage
                    if (VERBOSE>=1) {fprintf(stderr, "  (larger coverage too low.)\n");}
                    continue;
                }
                if (cov2<20){  // early termination, no significant coverage diff
                    if ((float)cov2/cov1<2 && (cov2-cov1<8)){
                        if (VERBOSE>=1) {fprintf(stderr, "  (diff insignificant, type 1. (%d vs %d))\n", (int)cov1, (int) cov2);}
                        continue;
                    }
                }else{
                    if ((float)cov2/cov1<2){
                        if (VERBOSE>=1) {fprintf(stderr, "  (diff insignificant, type 2. (%d vs %d))\n", (int)cov1, (int) cov2);}
                        continue;
                    }
                }

                ret = hamt_ugarc_detect_circle(ug, sg, vu<<1|dir, wu, v2vu_directed, sanmask);
                if (ret){
                    nb_cut++;
                    asg_arc_del(sg, v, w, 1);
                    asg_arc_del(sg, w^1, v^1, 1);
                    if (VERBOSE>=1){
                        fprintf(stderr, "[D::%s]dropped a link between utg%.6d and utg%.6d\n", __func__, (int)(vu)+1, (int)(wu>>1)+1);
                        fprintf(stderr, "        (coverage small: %d vs large: %d)\n", cov1, cov2);
                    }
                }else{
                    if (VERBOSE>=1) {fprintf(stderr, "  (circle not detected.)\n");}
                }
            }
        }
    }

    fprintf(stderr, "[M::%s] dropped %d links.\n", __func__, nb_cut);
    if (nb_cut){
        asg_cleanup(sg);
    }

    ma_ug_destroy(ug);
    free(v2vu_directed);
    free(sanmask);
    free(primary_flag);
}