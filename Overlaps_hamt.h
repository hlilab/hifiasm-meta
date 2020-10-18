#ifndef __OVERLAPS_HAMT__
#define __OVERLAPS_HAMT__
#include <stdint.h>
#include <assert.h>
#include "Overlaps.h"
#include "Process_Read.h"



void write_debug_assembly_graph(asg_t *sg, All_reads *rs, char* read_file_name);
int hamt_asgutil_detect_strict_single_long_path(asg_t *sg, uint32_t begNode, uint32_t *endNode, 
                                                int threshold_nodes, uint64_t threshold_length);
int hamt_asgutil_detect_dangling_circle(const asg_t *g, uint32_t v0, uint32_t w, uint8_t *buf_traversal);
int hamt_asgarc_util_get_the_one_target(asg_t *g, uint32_t v, int include_del_seq, int include_del_arc);
int hamt_asgarc_util_countSinglePath(asg_t *g, uint32_t v0, int include_del_seq, int include_del_arc, 
                                     int max_arcs, int max_length, asg64_v *a, int *traversal_dist);
int hamt_asgutil_calc_singlepath_cov(asg_t *g, uint32_t v);
void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length);
void hamt_asg_arc_del_by_readcov_circle_aware(asg_t *g);



// void hamt_asgarc_SCCcovCut(asg_t *g);
// void hamt_utg_SCCcovCut(asg_t *sg,
//                         const ma_sub_t* coverage_cut,   // for utg coverage
//                         ma_hit_t_alloc* sources, R_to_U* ruIndex);  // for utg coverage
// void hamt_usg_SCCcovCut(asg_t *sg0,
//                         const ma_sub_t* coverage_cut,   // for utg coverage
//                         ma_hit_t_alloc* sources, R_to_U* ruIndex);  // for utg coverage
void hamt_ugarc_covcut_danglingCircle(asg_t *sg,
                                        const ma_sub_t* coverage_cut,   // for utg coverage
                                        ma_hit_t_alloc* sources, R_to_U* ruIndex);  // for utg coverage
void hamt_asgarc_ugCovCutSCC(asg_t *sg,
                        const ma_sub_t* coverage_cut,   // for utg coverage
                        ma_hit_t_alloc* sources, R_to_U* ruIndex);  // for utg coverage

#endif // __OVERLAPS_HAMT__