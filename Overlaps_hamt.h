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

// interface: pre-contig-gen cleaning
void hamt_ug_prectgTopoClean(asg_t *sg, 
                            const ma_sub_t* coverage_cut, ma_hit_t_alloc* sources, R_to_U* ruIndex,
                            int base_label, int alt_label, int is_hard_drop);
void hamt_ug_prectg_rescueShortCircuit(asg_t *sg, ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, R_to_U* ruIndex,
                                        const ma_sub_t* coverage_cut, int base_label);
void hamt_ug_prectg_rescueShortCircuit_simpleAggressive(asg_t *sg, ma_ug_t *ug, 
                                                        ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources,
                                                        const ma_sub_t* coverage_cut, 
                                                        int base_label);
void hamt_ug_prectg_rescueLongUtg(asg_t *sg, 
                                    ma_hit_t_alloc *sources, ma_hit_t_alloc *reverse_sources, R_to_U* ruIndex,
                                    const ma_sub_t* coverage_cut);

int hamt_ug_prectg_resolve_complex_bubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop, int max_length);
int hamt_ug_resolve_oneMultiLeafSoapBubble(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);

// exp
int hamt_ug_resolveTangles(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label);

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










// tmp
int hamt_ug_pop_bubble(asg_t *sg,ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_pop_miscbubble(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_pop_miscbubble_aggressive(asg_t *sg, ma_ug_t *ug, int base_label);
int hamt_ug_pop_terminalSmallTip(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
int hamt_ug_pop_tinyUnevenCircle(asg_t *sg, ma_ug_t *ug, int base_label, int alt_label, int is_hard_drop);
#endif // __OVERLAPS_HAMT__