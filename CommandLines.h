#ifndef __COMMAND_LINE_PARSER__
#define __COMMAND_LINE_PARSER__

#include <pthread.h>

#define HA_VERSION "0.13-r308"
#define HAMT_VERSION "0.3-r063.2"


// #define VERBOSE 1
extern int VERBOSE;  // expose to cli

#define HA_F_NO_HPC          0x1
#define HA_F_NO_KMER_FLT     0x2
#define HA_F_VERBOSE_GFA     0x4
#define HA_F_WRITE_EC        0x8
#define HA_F_WRITE_PAF       0x10
#define HA_F_SKIP_TRIOBIN    0x20
#define HA_F_PURGE_CONTAIN   0x40
#define HA_F_PURGE_JOIN      0x80
#define HA_F_BAN_POST_JOIN   0x100
#define HA_F_BAN_ASSEMBLY    0x200
#define HA_F_HIGH_HET        0x400

#define HA_MIN_OV_DIFF       0.02 // min sequence divergence in an overlap


typedef struct {
	int flag;
    int num_reads;
    char** read_file_names;
    char* output_file_name;
    char* required_read_name;
	char *fn_bin_yak[2];
	char *fn_bin_list[2];
	char *extract_list;
	int extract_iter;
    int thread_num;
    int k_mer_length;
	int mz_win;
	int bf_shift;
	double high_factor; // coverage cutoff set to high_factor*hom_cov
	double max_ov_diff_ec;
	double max_ov_diff_final;
	int hom_cov;
    int het_cov;
	int max_n_chain; // fall-back max number of chains to consider
	int min_hist_kmer_cnt;
    int load_index_from_disk;
    int write_index_to_disk;
    int number_of_round;
    int adapterLen;
    int clean_round;
    int roundID;
    int max_hang_Len;
    int gap_fuzz;
    int min_overlap_Len;
    int min_overlap_coverage;
    int max_short_tip;
    int min_cnt;
    int mid_cnt;
    int purge_level_primary;
    int purge_level_trio;
    int purge_overlap_len;
    int recover_atg_cov_min;
    int recover_atg_cov_max;
    int hom_global_coverage;
    int bed_inconsist_rate;

    float max_hang_rate;
    float min_drop_rate;
    float max_drop_rate;
    float purge_simi_rate;

    long long small_pop_bubble_size;
    long long large_pop_bubble_size;
    long long num_bases;
    long long num_corrected_bases;
    long long num_recorrected_bases;
	long long mem_buf;
    long long coverage;

    // hamt
    int is_disable_read_selection;
    int mode_read_kmer_profile;
    int readselection_sort_order;  // experimental, 1 for smallestFirst, 2 for largestFirst, 0 to force disable it (note that we still go through loading all reads + sorting, just don't use the info when annotation mask_readnorm)
    char *bin_base_name;
    int is_ignore_ovlp_cnt;  // experimental, do preovec read selection even if the whole read set looks practical
    int preovec_coverage;
    int is_dump_read_mask;
    int is_dump_read_names;
    int is_use_exp_graph_cleaning;
    int is_dump_ovec_error_count;
    int lowq_thre_10;
    int lowq_thre_5;
    int lowq_thre_3;
    int write_debug_gfa;
    int write_new_graph_bins;
    int use_ha_bin;  // 

    int is_final_round;
    int is_mode_low_cov;
    int is_dump_relevant_reads; FILE *fp_relevant_reads;
    // hifiasm_argcv_t *argcv;

    // hamt, graph cleaning control
    int gc_superbubble_tig_max_length;
    int gc_tangle_max_tig;
    int mode_coasm;
    int is_aggressive;
    int do_probe_gfa;
    // hamt, multiprocessing context for get_specific_overlap's binary search 
    int get_specific_overlap_is_use_bf;

    // end of hamt
    
} hifiasm_opt_t;

typedef struct {
    int ha_argc;
    char **ha_argv;
} hifiasm_argcv_t;  // hamt

extern hifiasm_opt_t asm_opt;
extern hifiasm_argcv_t asm_argcv;

void init_opt(hifiasm_opt_t* asm_opt);
void destory_opt(hifiasm_opt_t* asm_opt);
void ha_opt_reset_to_round(hifiasm_opt_t* asm_opt, int round);
void ha_opt_update_cov(hifiasm_opt_t *opt, int hom_cov);
int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt);
double Get_T(void);

static inline int ha_opt_triobin(const hifiasm_opt_t *opt)
{
	return ((opt->fn_bin_yak[0] && opt->fn_bin_yak[1]) || (opt->fn_bin_list[0] && opt->fn_bin_list[1]));
}

#endif
