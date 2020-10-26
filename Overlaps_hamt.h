#ifndef __OVERLAPS_HAMT__
#define __OVERLAPS_HAMT__
#include <stdint.h>
#include <assert.h>
#include "Overlaps.h"
#include "Process_Read.h"

// debug / helper
void write_debug_assembly_graph(asg_t *sg, All_reads *rs, char* read_file_name);
void hamt_collect_utg_coverage(asg_t *sg, ma_ug_t *ug, 
                                const ma_sub_t* coverage_cut,
                                ma_hit_t_alloc* sources, R_to_U* ruIndex);
void hamt_destroy_utg_coverage(ma_ug_t *ug);
void hamt_ug_regen(asg_t *sg, ma_ug_t **ug,
                    const ma_sub_t* coverage_cut,
                    ma_hit_t_alloc* sources, R_to_U* ruIndex);


// exp routines 
void hamt_asgarc_drop_tips_and_bubbles(ma_hit_t_alloc* sources, asg_t *g, int max_arcs, int max_length);
void hamt_ug_covCutByBridges(asg_t *sg, ma_ug_t *ug);
void hamt_asgarc_ugCovCutDFSCircle(asg_t *sg, ma_ug_t *ug);
void hamt_asgarc_ugTreatMultiLeaf(asg_t *sg, ma_ug_t *ug, int threshold_l);

// interface
void hamt_circle_cleaning(asg_t *sg, ma_ug_t *ug);
void hamt_clean_shared_seq(asg_t *sg, ma_ug_t *ug);
void hamt_ug_pop_bubble(asg_t *sg, ma_ug_t *ug);
void hamt_ug_pop_miscbubble(asg_t *sg, ma_ug_t *ug);
void hamt_ug_pop_simpleInvertBubble(asg_t *sg, ma_ug_t *ug);


#endif // __OVERLAPS_HAMT__