#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <sys/time.h>
#include "CommandLines.h"
#include "ketopt.h"
#include "gitcommit.h"  // if gone, use stamp

#define DEFAULT_OUTPUT "hifiasm_meta.asm"

hifiasm_opt_t asm_opt;
hifiasm_argcv_t asm_argcv;  // hamt
int VERBOSE = 0; // expose to cli

static ko_longopt_t long_options[] = {
	{ "version",       ko_no_argument, 300 },
	{ "dbg-gfa",       ko_no_argument, 301 },  // write assembly graph bin files
	{ "write-paf",     ko_no_argument, 302 },
	{ "write-ec",      ko_no_argument, 303 },  // write error-corrected reads
	{ "skip-triobin",  ko_no_argument, 304 },
	{ "max-od-ec",     ko_required_argument, 305 },
	{ "max-od-final",  ko_required_argument, 306 },
	{ "ex-list",       ko_required_argument, 307 },
	{ "ex-iter",       ko_required_argument, 308 },
    { "purge-cov",     ko_required_argument, 309 },
    { "pri-range",     ko_required_argument, 310 },
    { "high-het",      ko_no_argument, 311 },

    // hamt debug/probing modules
    { "read-kmer-profile", ko_no_argument, 400},  // write per-read kmer frequency profiles
    { "dump-ovec-cnt", ko_no_argument, 402},  // dump per read number of corrected bases
    { "dump-read-mask", ko_no_argument, 406},  // dump the read selection mask from bin files; only effetive with -B (hamt)
    { "dump-read-names", ko_no_argument, 407},  // dump names of the selected reads
    { "force-preovec", ko_no_argument, 408}, // ignore 1st heuristic (which could've kept all reads), do preovec read selection based on lowq given
	
    { "lowq-10", ko_required_argument, 409}, // lower 10% quantile threshold
    { "lowq-5", ko_required_argument, 410}, // lower 5% quantile threshold
    { "lowq-3", ko_required_argument, 411}, // lower 3% quantile threshold

    { "inter-gfa", ko_no_argument, 413},  // write intermediate gfa files
    { "ban-meta", ko_no_argument, 414},  // use stable hifiasm route
    { "lowcov", ko_no_argument, 415},  // input has very low coverage, copy het overlaps to hom (experimental)
    { "refresh-dbg-gfa", ko_no_argument, 416},
    { "dump-all-ovlp", ko_no_argument, 417},

    { "gc-sb-max", ko_no_argument, 418},  // graph cleaning, max unitig length in superbubbles
    { "force-rs", ko_no_argument, 419},  // aka force-preovec, better named
    { "probe-gfa", ko_optional_argument, 420},  // compile-time probe of debug gfa, can pass prefix. MUST use = (gnu getopt behavior)
    { "use-ha-bin", ko_no_argument, 421}, // use hifiasm bin files - will ignore hamt-specific files and use placeholders.
    { "noch", ko_no_argument, 422},  // disable contained reads sparing heuristics
    { "write-binning", ko_no_argument, 423 },  // writes binning fasta files, in addition to the tsv

    { "tsne-perp", ko_required_argument, 424 },  // perplexity of tsne used in binning
    { "tsne-seed", ko_required_argument, 425 },  // 
    { "tsne-neighdist", ko_required_argument, 426 },  // 
    // end of hamt

    { "lowQ",          ko_required_argument, 312 },
	{ "min-hist-cnt",  ko_required_argument, 313 },
	{ 0, 0, 0 }
};

double Get_T(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

void Print_H(hifiasm_opt_t* asm_opt)
{
    fprintf(stderr, "Usage: hifiasm_meta %s (hifiasm code base %s)\n", HAMT_VERSION, HA_VERSION);
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Input/Output:\n");
    fprintf(stderr, "    -o STR      prefix of output files [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -B STR      prefix of bin files, if it's different from -o [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -i          ignore saved read correction and overlaps\n");
    fprintf(stderr, "    -t INT      number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    -z INT      length of adapters that should be removed [%d]\n", asm_opt->adapterLen);
    fprintf(stderr, "    --version   show version number\n");
    fprintf(stderr, "  Read selection:\n");
    // fprintf(stderr, "    -S          enable read selection.\n");
    fprintf(stderr, "    --force-rs\n");
    fprintf(stderr, "                enable and force read selection.\n");
    fprintf(stderr, "    --lowq-10\n");
    fprintf(stderr, "                lower 10%% runtime kmer frequency threshold. [%d]\n", asm_opt->lowq_thre_10);
    fprintf(stderr, "    --lowq-5\n");
    fprintf(stderr, "                lower 5%% runtime kmer frequency threshold. [%d]\n", asm_opt->lowq_thre_5);
    fprintf(stderr, "    --lowq-3\n");
    fprintf(stderr, "                lower 3%% runtime kmer frequency threshold. [%d]\n", asm_opt->lowq_thre_3);
    
	fprintf(stderr, "  Overlap/Error correction:\n");
    fprintf(stderr, "    -k INT      k-mer length (must be <64) [%d]\n", asm_opt->k_mer_length);
	fprintf(stderr, "    -w INT      minimizer window size [%d]\n", asm_opt->mz_win);
	fprintf(stderr, "    -f INT      number of bits for bloom filter; 0 to disable [%d]\n", asm_opt->bf_shift);
	fprintf(stderr, "    -D FLOAT    drop k-mers occurring >FLOAT*coverage times [%.1f]\n", asm_opt->high_factor);
	fprintf(stderr, "    -N INT      consider up to max(-D*coverage,-N) overlaps for each oriented read [%d]\n", asm_opt->max_n_chain);
    fprintf(stderr, "    -r INT      round of correction [%d]\n", asm_opt->number_of_round);
    fprintf(stderr, "  Assembly:\n");
    fprintf(stderr, "    -a INT      round of assembly cleaning [%d]\n", asm_opt->clean_round);
    fprintf(stderr, "    -n INT      remove tip unitigs composed of <=INT reads [%d]\n", asm_opt->max_short_tip);
    fprintf(stderr, "    -x FLOAT    max overlap drop ratio [%.2g]\n", asm_opt->max_drop_rate);
    fprintf(stderr, "    -y FLOAT    min overlap drop ratio [%.2g]\n", asm_opt->min_drop_rate);
    fprintf(stderr, "    --noch      disable contained reads sparing heuristics.\n");
    fprintf(stderr, "  Auxiliary:\n");
    fprintf(stderr, "    -e          ban assembly, i.e. terminate before generating string graph\n");
    fprintf(stderr, "    --write-paf dump overlaps (paf).\n");
    fprintf(stderr, "    --dump-all-ovlp\n");
    fprintf(stderr, "                dump all overlaps ever calculated in the final overlapping (paf).\n");
    fprintf(stderr, "    --write-ec\n");
    fprintf(stderr, "                dump error corrected reads (fasta).\n");


    fprintf(stderr, "Example: ./hifiasm_meta -o asm -t 32 asm.fq.gz 2>log\n");
    fprintf(stderr, "See `man ./hifiasm_meta.1` for detailed descriptions command-line options.\n");
}

void init_opt(hifiasm_opt_t* asm_opt)
{
	memset(asm_opt, 0, sizeof(hifiasm_opt_t));
	asm_opt->flag = 0;
    asm_opt->coverage = -1;
    asm_opt->num_reads = 0;
    asm_opt->read_file_names = NULL;
    asm_opt->output_file_name = (char*)(DEFAULT_OUTPUT);
    asm_opt->required_read_name = NULL;
    asm_opt->thread_num = 1;
    asm_opt->k_mer_length = 51;
	asm_opt->mz_win = 51;
	asm_opt->bf_shift = 37;
	asm_opt->high_factor = 5.0;
	asm_opt->max_ov_diff_ec = 0.04;
	asm_opt->max_ov_diff_final = 0.03;
	asm_opt->hom_cov = 20;
    asm_opt->het_cov = -1024;
	asm_opt->max_n_chain = 100;
	asm_opt->min_hist_kmer_cnt = 5;
    asm_opt->load_index_from_disk = 1;
    asm_opt->write_index_to_disk = 1;
    asm_opt->number_of_round = 3;
    asm_opt->adapterLen = 0;
    asm_opt->clean_round = 4;
    asm_opt->small_pop_bubble_size = 100000;
    asm_opt->large_pop_bubble_size = 10000000;
    asm_opt->min_drop_rate = 0.2;
    asm_opt->max_drop_rate = 0.8;
    asm_opt->max_hang_Len = 1000;
    asm_opt->max_hang_rate = 0.8;
    asm_opt->gap_fuzz = 1000;
    asm_opt->min_overlap_Len = 50;
    asm_opt->min_overlap_coverage = 0;
    asm_opt->max_short_tip = 3;
    asm_opt->min_cnt = 2;
    asm_opt->mid_cnt = 5;
    asm_opt->purge_level_primary = 2;
    asm_opt->purge_level_trio = 0;
    asm_opt->purge_simi_rate = 0.75;
    asm_opt->purge_overlap_len = 1;
    asm_opt->recover_atg_cov_min = -1024;
    asm_opt->recover_atg_cov_max = INT_MAX;
    asm_opt->hom_global_coverage = -1;

    // hamt
    asm_opt->is_disable_read_selection = 1;
    asm_opt->mode_read_kmer_profile = 0;
    asm_opt->readselection_sort_order = 1;  // smallest first
    asm_opt->bin_base_name = (char*)(DEFAULT_OUTPUT);
    asm_opt->probebin_base_name = 0;
    asm_opt->preovec_coverage = 150;
    asm_opt->is_ignore_ovlp_cnt = 0;
    asm_opt->is_dump_read_mask = 0;
    asm_opt->is_dump_read_names = 0;
    asm_opt->is_use_exp_graph_cleaning = 1;
    asm_opt->is_dump_ovec_error_count = 0;
    asm_opt->lowq_thre_10 = 50;  // lower 10% quantile runtime kmer frequency
    asm_opt->lowq_thre_5 = 50;  // lower 5% quantile runtime kmer frequency
    asm_opt->lowq_thre_3 = 10;  // lower 5% quantile runtime kmer frequency
    asm_opt->write_debug_gfa = 0;  // disable
    asm_opt->is_dump_relevant_reads = 0;
    asm_opt->fp_relevant_reads = NULL;
    asm_opt->is_mode_low_cov = 0;
    asm_opt->write_new_graph_bins = 0;
    asm_opt->gc_superbubble_tig_max_length = 100000;
    asm_opt->gc_tangle_max_tig = 500;  // set to -1 to disable the limit
    asm_opt->is_aggressive = 0; 
    asm_opt->use_ha_bin = 0;
    asm_opt->no_containedreads_heuristics = 0;
    asm_opt->write_binning_fasta = 0;  // defaults to only write a tsv for binning
    asm_opt->tsne_perplexity = 50;
    asm_opt->tsne_randomseed = 42; 
    asm_opt->tsne_neigh_dist = 0.25;
    // end of hamt
    asm_opt->bed_inconsist_rate = 0;  // hamt: disable
}

void destory_opt(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names != NULL)
    {
        free(asm_opt->read_file_names);
    }
}

void ha_opt_reset_to_round(hifiasm_opt_t* asm_opt, int round)
{
    asm_opt->num_bases = 0;
    asm_opt->num_corrected_bases = 0;
    asm_opt->num_recorrected_bases = 0;
	asm_opt->mem_buf = 0;
    asm_opt->roundID = round;
}

void ha_opt_update_cov(hifiasm_opt_t *opt, int hom_cov)
{
	int max_n_chain = (int)(hom_cov * opt->high_factor + .499);
	opt->hom_cov = hom_cov;
	if (opt->max_n_chain < max_n_chain)
		opt->max_n_chain = max_n_chain;
	fprintf(stderr, "[M::%s] updated max_n_chain to %d\n", __func__, opt->max_n_chain);
}

static int check_file(char* name, const char* opt)
{
    if(!name)
    {
        fprintf(stderr, "[ERROR] file does not exist (-%s)\n", opt);
        return 0;
    } 
    FILE* is_exist = NULL;
    is_exist = fopen(name,"r");
    if(!is_exist)
    {
        fprintf(stderr, "[ERROR] %s does not exist (-%s)\n", name, opt);
        return 0;
    } 

    fclose(is_exist);
    return 1;
}

int check_option(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names == NULL || asm_opt->num_reads == 0)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if(asm_opt->output_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name (-o)\n");
        return 0;
    }

    if(asm_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0 (-t)\n");
        return 0;
    }


    if(asm_opt->number_of_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for correction must be > 0 (-r)\n");
        return 0;
    }

    if(asm_opt->clean_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for assembly cleaning must be > 0 (-a)\n");
        return 0;
    }

    if(asm_opt->adapterLen < 0)
    {
        fprintf(stderr, "[ERROR] the length of removed adapters must be >= 0 (-z)\n");
        return 0;
    }


    if(asm_opt->k_mer_length >= 64)
    {
        fprintf(stderr, "[ERROR] the length of k_mer must be < 64 (-k)\n");
        return 0;
    }


    if(asm_opt->max_drop_rate < 0 || asm_opt->max_drop_rate >= 1 ) 
    {
        fprintf(stderr, "[ERROR] max overlap drop ratio must be [0.0, 1.0) (-x)\n");
        return 0;
    }

    
    if(asm_opt->min_drop_rate < 0 || asm_opt->min_drop_rate >= 1)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be [0.0, 1.0) (-y)\n");
        return 0;
    }

    if(asm_opt->max_drop_rate <= asm_opt->min_drop_rate)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be less than max overlap drop ratio (-x/-y)\n");
        return 0;
    }

    if(asm_opt->small_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped small bubbles must be >= 0 (-p)\n");
        return 0;
    }

    if(asm_opt->large_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped large bubbles must be >= 0 (-m)\n");
        return 0;
    }

    if(asm_opt->max_hang_Len < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->max_hang_rate < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_rate must be >= 0\n");
        return 0;
    }

    if(asm_opt->gap_fuzz < 0)
    {
        fprintf(stderr, "[ERROR] gap_fuzz must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_Len < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_coverage < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_coverage must be >= 0\n");
        return 0;
    }

	if (asm_opt->max_ov_diff_ec < asm_opt->max_ov_diff_final) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than max_ov_diff_final\n");
		return 0;
	}

	if (asm_opt->max_ov_diff_ec < HA_MIN_OV_DIFF) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than %g\n", HA_MIN_OV_DIFF);
		return 0;
	}

    if(asm_opt->max_short_tip < 0)
    {
        fprintf(stderr, "[ERROR] the length of removal tips must be >= 0 (-n)\n");
        return 0;
    }

    if(asm_opt->purge_level_primary < 0 || asm_opt->purge_level_primary > 2)
    {
        fprintf(stderr, "[ERROR] the level of purge-dup should be [0, 2] (-l)\n");
        return 0;
    }

    if(ha_opt_triobin(asm_opt) && ((asm_opt->purge_level_trio < 0 || asm_opt->purge_level_trio > 1)))
    {
        fprintf(stderr, "[ERROR] the level of purge-dup for trio should be [0, 1] (-l)\n");
        return 0;
    }

    if(asm_opt->hom_global_coverage < 0 && asm_opt->hom_global_coverage != -1)
    {
        fprintf(stderr, "[ERROR] purge duplication coverage threshold should be >= 0 (--purge-cov)\n");
        return 0;
    }

    if(asm_opt->bed_inconsist_rate < 0 || asm_opt->bed_inconsist_rate > 100)
    {
        fprintf(stderr, "[ERROR] inconsistency rate should be [0, 100] (--pb-range)\n");
        return 0;
    }


    if(asm_opt->fn_bin_yak[0] != NULL && check_file(asm_opt->fn_bin_yak[0], "YAK1") == 0) return 0;
    if(asm_opt->fn_bin_yak[1] != NULL && check_file(asm_opt->fn_bin_yak[1], "YAK2") == 0) return 0;
    if(asm_opt->fn_bin_list[0] != NULL && check_file(asm_opt->fn_bin_list[0], "LIST1") == 0) return 0;
    if(asm_opt->fn_bin_list[1] != NULL && check_file(asm_opt->fn_bin_list[1], "LIST2") == 0) return 0;
    if(asm_opt->required_read_name != NULL && check_file(asm_opt->required_read_name, "b") == 0) return 0;
    // fprintf(stderr, "input file num: %d\n", asm_opt->num_reads);
    // fprintf(stderr, "output file: %s\n", asm_opt->output_file_name);
    // fprintf(stderr, "number of threads: %d\n", asm_opt->thread_num);
    // fprintf(stderr, "number of rounds for correction: %d\n", asm_opt->number_of_round);
    // fprintf(stderr, "number of rounds for assembly cleaning: %d\n", asm_opt->clean_round);
    // fprintf(stderr, "length of removed adapters: %d\n", asm_opt->adapterLen);
    // fprintf(stderr, "length of k_mer: %d\n", asm_opt->k_mer_length);
    // fprintf(stderr, "min overlap drop ratio: %.2g\n", asm_opt->min_drop_rate);
    // fprintf(stderr, "max overlap drop ratio: %.2g\n", asm_opt->max_drop_rate);
    // fprintf(stderr, "size of popped small bubbles: %lld\n", asm_opt->small_pop_bubble_size);
    // fprintf(stderr, "size of popped large bubbles: %lld\n", asm_opt->large_pop_bubble_size);
    // fprintf(stderr, "small removed unitig threshold: %d\n", asm_opt->max_short_tip);
    // fprintf(stderr, "small removed unitig threshold: %d\n", asm_opt->max_short_tip);
    // fprintf(stderr, "min_cnt: %d\n", asm_opt->min_cnt);
    // fprintf(stderr, "mid_cnt: %d\n", asm_opt->mid_cnt);
    // fprintf(stderr, "purge_level_primary: %d\n", asm_opt->purge_level_primary);
    // fprintf(stderr, "purge_level_trio: %d\n", asm_opt->purge_level_trio);
    // fprintf(stderr, "purge_simi_rate: %f\n", asm_opt->purge_simi_rate);
    // fprintf(stderr, "purge_overlap_len: %d\n", asm_opt->purge_overlap_len);

    ///// hamt /////
    // if (asm_opt->lowq_thre_10<=1) {fprintf(stderr, "[E::%s] lowq-10 threshold too small.\n", __func__); return 0;}

    return 1;
}

void get_queries(int argc, char *argv[], ketopt_t* opt, hifiasm_opt_t* asm_opt)
{
    if(opt->ind == argc)
    {
        return;
    }

    asm_opt->num_reads = argc - opt->ind;
    asm_opt->read_file_names = (char**)malloc(sizeof(char*)*asm_opt->num_reads);
    
    long long i;
    gzFile dfp;
    for (i = 0; i < asm_opt->num_reads; i++)
    {
        asm_opt->read_file_names[i] = argv[i + opt->ind];
        dfp = gzopen(asm_opt->read_file_names[i], "r");
        if (dfp == 0)
        {
            fprintf(stderr, "[ERROR] Cannot find the input read file: %s\n", 
                    asm_opt->read_file_names[i]);
		    exit(0);
        }
        gzclose(dfp);
    }
}


int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    asm_argcv.ha_argc = argc;
    asm_argcv.ha_argv = argv;    
    // asm_opt->argcv = &asm_argcv;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o:k:w:m:n:r:a:b:z:x:y:p:c:d:M:P:if:D:FN:1:2:3:4:l:s:O:eu:VSB:A", long_options)) >= 0) {
        if (c == 'h')
        {
            Print_H(asm_opt);
            return 0;
        } 
        else if (c == 'v' || c == 300)
        {
            fprintf(stderr, "ha base version: %s\n", HA_VERSION);
            fprintf(stderr, "hamt version: %s\n", HAMT_VERSION);
            return 0;
        }
        // vanilla: abcdef hi klmnopqrstuvwxyz
        // vanilla:    D F      MNOP               
        // meta   :
        // meta   : AB                S  V 
		else if (c == 'f') asm_opt->bf_shift = atoi(opt.arg);
        else if (c == 't') asm_opt->thread_num = atoi(opt.arg); 
        else if (c == 'o') {
            asm_opt->output_file_name = opt.arg;
            if (strcmp(asm_opt->bin_base_name, DEFAULT_OUTPUT)==0) asm_opt->bin_base_name = opt.arg;
        }
        else if (c == 'r') asm_opt->number_of_round = atoi(opt.arg);
        else if (c == 'k') asm_opt->k_mer_length = atoi(opt.arg);
        else if (c == 'i') asm_opt->load_index_from_disk = 0; 
        else if (c == 'w') asm_opt->mz_win = atoi(opt.arg);
		else if (c == 'D') asm_opt->high_factor = atof(opt.arg);
		else if (c == 'F') asm_opt->flag |= HA_F_NO_KMER_FLT;
		else if (c == 'N') asm_opt->max_n_chain = atoi(opt.arg);
        else if (c == 'a') asm_opt->clean_round = atoi(opt.arg); 
        else if (c == 'z') asm_opt->adapterLen = atoi(opt.arg);
        else if (c == 'b') asm_opt->required_read_name = opt.arg;
        else if (c == 'c') asm_opt->min_cnt = atoi(opt.arg);
        else if (c == 'd') asm_opt->mid_cnt = atoi(opt.arg);
        else if (c == '1' || c == 'P') asm_opt->fn_bin_yak[0] = opt.arg; // -P/-M reserved for backward compatibility
        else if (c == '2' || c == 'M') asm_opt->fn_bin_yak[1] = opt.arg;
        else if (c == '3') asm_opt->fn_bin_list[0] = opt.arg;
        else if (c == '4') asm_opt->fn_bin_list[1] = opt.arg;
        else if (c == 'x') asm_opt->max_drop_rate = atof(opt.arg);
        else if (c == 'y') asm_opt->min_drop_rate = atof(opt.arg);
        else if (c == 'p') asm_opt->small_pop_bubble_size = atoll(opt.arg);
        else if (c == 'm') asm_opt->large_pop_bubble_size = atoll(opt.arg);
        else if (c == 'n') asm_opt->max_short_tip = atoll(opt.arg);
        else if (c == 'e') asm_opt->flag |= HA_F_BAN_ASSEMBLY;
        else if (c == 'u') asm_opt->flag |= HA_F_BAN_POST_JOIN;
        // hamt
        else if (c == 'A') {
            asm_opt->is_aggressive = 1;
            asm_opt->gc_tangle_max_tig = 1000;
        }
        else if (c == 'V') VERBOSE += 1;  // 1 will print out ha's debug and a few others, 1+ will print ovlp read skip info for each read
        else if (c == 'S') {
            asm_opt->is_disable_read_selection = 0; 
            fprintf(stderr, "[M::%s] Read selection enabled.\n", __func__);
        }
        else if (c == 'B') {asm_opt->bin_base_name = opt.arg; fprintf(stderr, "[M::%s] Use bin files under the name %s\n", __func__, opt.arg);}  // using bin files from another location and/or under different name
        else if (c == 400) {asm_opt->mode_read_kmer_profile = 1; fprintf(stderr, "DEBUG DUMP: get kmer frequency profile for every read.\n");} 
        else if (c == 402) {asm_opt->is_dump_ovec_error_count = 1; fprintf(stderr, "DEBUG DUMP: get ovec error counts\n");}
        else if (c == 406) {asm_opt->is_dump_read_mask = 1; fprintf(stderr, "DEBUG DUMP: will write read selection mask to file.\n");}
        else if (c == 407) {asm_opt->is_dump_read_names = 1; fprintf(stderr, "DEBUG DUMP: read names\n");}
        else if (c == 408 || c==419) {
            fprintf(stderr, "[M::%s] Forced pre-ovec read selection. Ignoring count of ovlp.\n", __func__);
            asm_opt->is_ignore_ovlp_cnt = 1;
            asm_opt->is_disable_read_selection = 0;
        }
        else if (c == 409) {
            fprintf(stderr, "[M::%s] Set lowq(10%%); note that without --force-preovec, this threshold have no effect if total number of overlaps is considered to be acceptable.\n", __func__);
            asm_opt->lowq_thre_10 = atoi(opt.arg);
            asm_opt->is_disable_read_selection = 0;
        }
        else if (c == 410) {
            fprintf(stderr, "[M::%s] Set lowq(5%%); note that without --force-preovec, this threshold have no effect if total number of overlaps is considered to be acceptable.\n", __func__);
            asm_opt->lowq_thre_5 = atoi(opt.arg);
            asm_opt->is_disable_read_selection = 0;
        }
        else if (c == 411) {
            fprintf(stderr, "[M::%s] Set lowq(3%%); note that without --force-preovec, this threshold have no effect if total number of overlaps is considered to be acceptable.\n", __func__);
            asm_opt->lowq_thre_3 = atoi(opt.arg);
            asm_opt->is_disable_read_selection = 0;
        }
        else if (c == 413) {
            fprintf(stderr, "[M::%s] will write intermediate gfa files; also set VERBOSE to 1\n", __func__);
            asm_opt->write_debug_gfa = 1;
            VERBOSE = 1;
        }
        else if (c == 414) {asm_opt->is_use_exp_graph_cleaning = 0;}
        else if (c == 415) {asm_opt->is_mode_low_cov = 1;}
        else if (c == 416) {asm_opt->write_new_graph_bins = 1;}
        else if (c == 417) {
            fprintf(stderr, "FILE DUMP: will dump all overlaps that's been considered before -N threshold is applied\n");
            fprintf(stderr, "     hint: use --write-ec to dump error-corrected reads.\n");
            fprintf(stderr, "           specify bin files from other dir with -B switch.\n");
            asm_opt->is_dump_relevant_reads = 1;
        }
        else if (c == 418) {asm_opt->gc_superbubble_tig_max_length = atoi(opt.arg);}
        else if (c == 420) {
            asm_opt->do_probe_gfa = 1;  
            if (opt.arg) {
                asm_opt->probebin_base_name = opt.arg;
                fprintf(stderr, "[M::%s] may use probe gfa bins at %s\n", __func__, opt.arg);
            }else{
                fprintf(stderr, "[M::%s] probe-gfa suffix will fall back to -o or -B, note that to pass a string you need the equal sign.\n", __func__);
            }
        }
        else if (c == 421) {asm_opt->use_ha_bin = 1; fprintf(stderr, "[M::%s] use hifiasm bin files\n", __func__);}
        else if (c == 422) {asm_opt->no_containedreads_heuristics = 1; 
                            fprintf(stderr, "[M::%s] contained reads sparing heuristics disabled.\n", __func__);}
        else if (c == 423) {asm_opt->write_binning_fasta = 1;}
        else if (c == 424) {asm_opt->tsne_perplexity= atoi(opt.arg);}
        else if (c == 425) {asm_opt->tsne_randomseed= atoi(opt.arg);}
        else if (c == 426) {asm_opt->tsne_neigh_dist= atof(opt.arg);}

        // end of hamt
		else if (c == 301) asm_opt->flag |= HA_F_VERBOSE_GFA;
		else if (c == 302) asm_opt->flag |= HA_F_WRITE_PAF;
		else if (c == 303) asm_opt->flag |= HA_F_WRITE_EC;
		else if (c == 304) asm_opt->flag |= HA_F_SKIP_TRIOBIN;
		else if (c == 305) asm_opt->max_ov_diff_ec = atof(opt.arg);
		else if (c == 306) asm_opt->max_ov_diff_final = atof(opt.arg);
		else if (c == 307) asm_opt->extract_list = opt.arg;
		else if (c == 308) asm_opt->extract_iter = atoi(opt.arg);
        else if (c == 309) asm_opt->hom_global_coverage = atoi(opt.arg);
        else if (c == 310)
        {
            char* s = NULL;
            asm_opt->recover_atg_cov_min = strtol(opt.arg, &s, 10);
			if (*s == ',') asm_opt->recover_atg_cov_max = strtol(s + 1, &s, 10);
            if(asm_opt->recover_atg_cov_min == -1 || asm_opt->recover_atg_cov_max == -1)
            {
                asm_opt->recover_atg_cov_min = asm_opt->recover_atg_cov_max = -1;
            }
        }
        else if (c == 311) asm_opt->flag |= HA_F_HIGH_HET;
        else if (c == 312) asm_opt->bed_inconsist_rate = atoi(opt.arg);
		else if (c == 313) asm_opt->min_hist_kmer_cnt = atoi(opt.arg);
        else if (c == 'l')
        {   ///0: disable purge_dup; 1: purge containment; 2: purge overlap
            asm_opt->purge_level_primary = asm_opt->purge_level_trio = atoi(opt.arg);
        }
        else if (c == 's') asm_opt->purge_simi_rate = atof(opt.arg);
        else if (c == 'O') asm_opt->purge_overlap_len = atoll(opt.arg);
        else if (c == ':') 
        {
			fprintf(stderr, "[ERROR] missing option argument in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		} 
        else if (c == '?') 
        {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		}
    }

    // swith order should not matter. need to process some cases here.
    if (asm_opt->do_probe_gfa){
        if (!asm_opt->probebin_base_name) {
            asm_opt->probebin_base_name = asm_opt->bin_base_name;
            fprintf(stderr, "[M::%s] update probe-gfa prefix to %s\n", __func__, 
                    asm_opt->bin_base_name);
        }
    }

    if (argc == opt.ind)
    {
        Print_H(asm_opt);
        return 0;
    }

    get_queries(argc, argv, &opt, asm_opt);

    return check_option(asm_opt);
}
