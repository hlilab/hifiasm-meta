#ifndef __OVERLAPS_HAMT__
#define __OVERLAPS_HAMT__
#include <stdint.h>
#include <assert.h>
#include "Overlaps.h"
#include "Process_Read.h"

typedef struct {
    char *a;
    int n, m;
} dbgmsg_t;
// for debug fprintf from kthreads
// (note: printf operates in multithreaded; might crash if race)
void hamt_dbgmsg_init(dbgmsg_t *h);
void hamt_dbgmsg_destroy(dbgmsg_t *h);
void hamt_dbgmsg_reset(dbgmsg_t *h);
void hamt_dbgmsg_append(dbgmsg_t *h, char *s, int l);
int hamt_dbgmsg_is_empty(dbgmsg_t *h);

// debug / helper
void write_debug_assembly_graph(asg_t *sg, All_reads *rs, char* read_file_name);
void hamt_collect_utg_coverage(asg_t *sg, ma_ug_t *ug, 
                                const ma_sub_t* coverage_cut,
                                ma_hit_t_alloc* sources, R_to_U* ruIndex);
void hamt_destroy_utg_coverage(ma_ug_t *ug);
ma_ug_t *hamt_ug_gen(asg_t *sg,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex, int flag);
void hamt_ug_destroy(ma_ug_t *ug);
void hamt_ug_regen(asg_t *sg, ma_ug_t **ug,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex, int label);
void hamtdebug_output_unitig_graph_ug(ma_ug_t *ug, char *base_fn, const char *suffix, int cleanID);
void hamt_ug_init_seq_vis(ma_ug_t *ug, uint8_t flag);
void hamt_asg_reset_seq_label(asg_t *sg, uint8_t flag);
void hamt_asg_reset_seq_vis(asg_t *sg, uint8_t flag);
void hamt_debug_dump(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources);

int hamt_debug_ug_random_cut_arcs(asg_t *sg, ma_ug_t *ug, int nb_cut);

// file output helper routines
int hamt_ug_util_BFS_markSubgraph(ma_ug_t *ug, int base_label);

// exp routines 
void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length);
void hamt_ug_covCutByBridges(asg_t *sg, ma_ug_t *ug, int base_label);
void hamt_asgarc_ugCovCutDFSCircle(asg_t *sg, ma_ug_t *ug, int base_label);
void hamt_asgarc_ugTreatMultiLeaf(asg_t *sg, ma_ug_t *ug, int threshold_l, int base_label, int alt_label, int is_hard_drop);


// interface: cleaning
int hamt_ugasg_cut_shortTips(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_cut_shortTips_arbitrary(asg_t *sg, ma_ug_t *ug, int max_length, int base_label);
int hamt_ug_covcut_falseVertexLoop(asg_t *sg, ma_ug_t *ug, int base_label);
void hamt_circle_cleaning(asg_t *sg, ma_ug_t *ug, int base_label);
void hamt_clean_shared_seq(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_basic_topoclean(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_basic_topoclean_simple(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_special_toploclean(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label);
int hamt_ug_pop_simpleInvertBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
void hamt_asgarc_ugCovCutDFSCircle_aggressive(asg_t *sg, ma_ug_t *ug, int base_label);

int hamt_ug_pop_simpleShortCut(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_oneutgCircleCut(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_oneutgCircleCut2(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_drop_midsizeTips(asg_t *sg, ma_ug_t *ug, int fold, int base_label);
int hamt_ug_drop_midsizeTips_aggressive(asg_t *sg, ma_ug_t *ug, float fold, int base_label);
int hamt_ug_resolve_small_multileaf_with_covcut(asg_t *sg, ma_ug_t *ug, int max_length, int fold, int base_label);

int hamt_ug_drop_transitive(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label);

int hamt_ug_drop_redundant_nodes(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label);
int hamt_ug_drop_redundant_nodes_bruteforce(asg_t *sg, ma_ug_t *ug, int size_limit_bp, int base_label, int verbose);

// interface: pre-contig-gen cleaning
void hamt_ug_prectgTopoClean(asg_t *sg, 
                            const ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex,
                            int base_label, int alt_label, int is_hard_drop);
int hamt_ug_prectg_rescueShortCircuit(asg_t *sg, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, R_to_U* ruIndex,
                                        const ma_sub_t* coverage_cut, int base_label);
int hamt_ug_prectg_rescueShortCircuit_simpleAggressive(asg_t *sg, ma_ug_t *ug, 
                                                        ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,
                                                        const ma_sub_t* coverage_cut, 
                                                        int base_label);
void hamt_ug_rescueLongUtg(asg_t *sg, 
                                    ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, R_to_U* ruIndex,
                                    const ma_sub_t* coverage_cut);
int hamt_ug_try_circularize(asg_t *sg, ma_ug_t *ug, 
                            ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,R_to_U* ruIndex,
                            const ma_sub_t* coverage_cut, int l_threshold);

int hamt_ug_prectg_resolve_complex_bubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop, int max_length);
int hamt_ug_resolve_oneMultiLeafSoapBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_pop_tinyFlatCircles(asg_t *sg, ma_ug_t *ug, int base_label);

// exp
int hamt_ug_resolveTangles(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label);
int hamt_ug_pop_unevenInvertBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label);

// (hap)
int hamt_ug_rescueLowCovHapGap_simple(asg_t *sg, ma_ug_t *ug, 
                              ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, 
                              const ma_sub_t* coverage_cut, uint64_t *readLen);
int hamt_ug_treatBifurcation_hapCovCut(asg_t *sg, ma_ug_t *ug, float covdiff_ratio, float haplo_ratio, 
                                        ma_hit_t_alloc *reverse_sources,
                                        int base_label, int alt_label);
void hamt_ughit_rescueLowCovHapGap(asg_t *sg, ma_ug_t *ug, 
                                         ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, const ma_sub_t* coverage_cut,
                                         long long n_read, uint64_t *readLen, int read_cov_threshold);
int hamt_ug_resolve_fake_haplotype_bifurcation(asg_t *sg, ma_ug_t *ug, int base_label,
                                               ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources);
int hamt_ug_resolve_fake_haplotype_bifurcation_aggressive(asg_t *sg, ma_ug_t *ug, int base_label,
                                               ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources);
void hamt_debug_get_diploid_info_about_all_branchings(ma_ug_t *ug, ma_hit_t_alloc *reverse_sources);
int hamt_ug_rescue_bifurTip(asg_t *sg, ma_ug_t *ug,int base_label,
                           ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, 
                           const ma_sub_t* coverage_cut);

// misc
int hamt_ug_cleanup_almost_circular(asg_t *sg, ma_ug_t *ug, int base_label);









// tmp
int hamt_ug_pop_bubble(asg_t *sg,ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_pop_miscbubble(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_pop_miscbubble_aggressive(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_pop_terminalSmallTip(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_pop_tinyUnevenCircle(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);

// coasm
int hamt_asg_arc_del_intersample_branching(asg_t *sg,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex);  // do not use
int hamt_ug_cut_very_short_multi_tip(asg_t *sg, ma_ug_t *ug, int nb_threshold);
int hamt_ug_drop_shorter_ovlp(asg_t *sg, ma_ug_t *ug, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources);


// 
int hamt_ug_finalprune(asg_t *sg, ma_ug_t *ug);
int hamt_ug_popLooseTangles(asg_t *sg, ma_ug_t *ug, int threshold_min_handle_length);
int hamt_ug_popLooseTangles_v2(asg_t *sg, ma_ug_t *ug, int max_step);


// binning experimental stuff

void hamt_dump_path_coverage_with_haplotype_info(ma_ug_t *ug, asg_t *read_g, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              R_to_U* ruIndex, const ma_sub_t* coverage_cut,
                              char *asm_prefix);
void hamt_dump_haplotig_pairs(ma_ug_t *ug, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              char *asm_prefix);
void hamt_update_coverage(ma_ug_t *ug, asg_t *read_g, 
                              ma_hit_t_alloc * sources, ma_hit_t_alloc *reverse_sources, 
                              R_to_U* ruIndex, const ma_sub_t* coverage_cut,
                              char *asm_prefix);





void hamt_utg_scc_testing(ma_ug_t *ug, int *labels);
void hamt_ug_get_all_elementary_circuits(ma_ug_t *ug);
void hamt_ug_opportunistic_elementary_circuits(asg_t *sg, ma_ug_t *ug, int n_thread);
int hamt_ug_delete_unconnected_single_read_contigs(asg_t *sg, ma_ug_t *ug);

#endif // __OVERLAPS_HAMT__