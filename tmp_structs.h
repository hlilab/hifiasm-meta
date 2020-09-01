#include <stdint.h>
#include "kseq.h"
#include "kvec.h"
#include "kthread.h"
#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, x|=(x)>>32, ++(x))
#define YAK_COUNTER_BITS 12
#define YAK_N_COUNTS (1<<YAK_COUNTER_BITS)

//#define CORRECT_THRESHOLD 0.70
#define CORRECT_THRESHOLD 0.60
///#define CORRECT_THRESHOLD_SECOND 0.55
#define CORRECT_THRESHOLD_HOMOPOLYMER 0.515
#define MIN_COVERAGE_THRESHOLD 3
#define CORRECT_INDEL_LENGTH 2
#define MISMATCH 1
#define INSERTION 2
#define DELETION 3

#define WINDOW_MAX_SIZE (WINDOW + (int)(1.0 / HA_MIN_OV_DIFF) + 3) // TODO: why 1/max_ov_diff?

//************************
//      Process_Read
//************************

typedef struct
{
    uint64_t x_id;
    uint64_t x_pos_s;
    uint64_t x_pos_e;
    uint8_t x_pos_strand;

    uint64_t y_id;
    uint64_t y_pos_s;
    uint64_t y_pos_e;
    uint8_t y_pos_strand;

	uint64_t matchLen;
	uint64_t totalLen;

} PAF;

typedef struct
{
    PAF* list;
    uint64_t size;
    uint64_t length;
} PAF_alloc;

typedef struct
{
    /**[0-1] bits are type:**/
    /**[2-31] bits are length**/
    uint32_t* record;
    uint32_t length;
	uint32_t size;

    char* lost_base;
    uint32_t lost_base_length;
	uint32_t lost_base_size;
	uint32_t new_length;
} Compressed_Cigar_record;

#define AMBIGU 0
#define FATHER 1
#define MOTHER 2
#define MIX_TRIO 3
#define NON_TRIO 4
#define DROP 5

typedef struct
{
	uint64_t** N_site;
	char* name;

	uint8_t** read_sperate;
	uint64_t* read_length; // alloc by `init_All_reads`
	uint64_t* read_size;  // alloc by `malloc_All_reads`
	uint8_t* trio_flag;

	uint64_t index_size;  // set by `init_All_reads`; hard-coded as 1000 by Process_Read.h

    ///name start pos in char* name
	uint64_t* name_index;  // alloc by `init_All_reads`
	uint64_t name_index_size;
	uint64_t total_reads;
	uint64_t total_reads_bases;
	uint64_t total_name_length;

	Compressed_Cigar_record* cigars; 
	Compressed_Cigar_record* second_round_cigar;

    ma_hit_t_alloc* paf;
    ma_hit_t_alloc* reverse_paf;
} All_reads;

typedef struct
{
	char* seq;
	long long length;
	long long size;
	long long RID;
} UC_Read;

//************************
//       Assembly
//************************
typedef struct {  // per thread data structure
	int is_final, save_ov;
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
	ha_abuf_t *ab;
	// error correction related buffers
	int64_t num_read_base, num_correct_base, num_recorrect_base;
	Cigar_record cigar1;
	Graph POA_Graph;
	Graph DAGCon;
	Correct_dumy correct;
	haplotype_evdience_alloc hap;
	Round2_alignment round2;
} ha_ovec_buf_t;

typedef struct {
	UC_Read g_read;
	int first_round_read_size;
	int second_round_read_size;
	char *first_round_read;
	char *second_round_read;
} ha_ecsave_buf_t;

//************************
//       htab
//************************

typedef struct {
	uint64_t x;  // 64 bit hash value of the <64bp minimizer/kmer
	uint64_t rid:28, pos:27, rev:1, span:8;
} ha_mz1_t;

typedef struct {
	uint64_t rid:28, pos:27, rev:1, span:8; // actually it is not necessary to keep span in the index
} ha_idxpos_t;

typedef struct {   // array of minimizers
    uint32_t n, m; 
    ha_mz1_t *a; 
} ha_mz1_v;  

struct ha_pt_s;
typedef struct ha_pt_s ha_pt_t;

struct ha_abuf_s;
typedef struct ha_abuf_s ha_abuf_t;

typedef struct {  // buffer for counting all kmers
// a linear buffer that stores the kmer hashes, ordering determined by how reads were loaded & parsed.
// used as a ensemble
	int n, m;
	uint64_t n_ins;
	uint64_t *a;  // the linear buffer for ct, each entry is a hash value of a kmer
	ha_mz1_t *b;  // the linear buffer for pt
} ch_buf_t;

typedef struct { // kmer counting; global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store;
	uint64_t n_seq;  // number of sequences; has to be less than 1<<28
	kseq_t *ks;  // will be there if reads aren't in mem
	UC_Read ucr;  // decompressed read
	ha_ct_t *ct;  // ensemble of minimizer count hashtables & bloom filters.
	ha_pt_t *pt;  // ensemble of minimizer position hashtables & ha_idxpost_t entries (uint64_t rid:28, pos:27, rev:1, span:8;)
	const All_reads *rs_in;
	All_reads *rs_out;  // if HAF_RS_WRITE_LEN or HAF_RS_WRITE_SEQ is set (not both!)
} pl_data_t;

typedef struct { // kmer counting; data structure for each step in kt_pipeline()
	pl_data_t *p;
	uint64_t n_seq0;  // rid of the first sequence (in the current state)?
	int n_seq, m_seq, sum_len, nk;  // sum_len: total loaded sequence length
	int *len;
	char **seq;  // sequences. not compressed.
	ha_mz1_v *mz_buf;  // per thread, array of minimizers
	ha_mz1_v *mz;  // per read holder
	ch_buf_t *buf;  // ensemble of linear kmer buffers; an array of ha_mz1_t entries w/ helpers. { int n, m; uint64_t n_ins; uint64_t *a; ha_mz1_t *b; }
} st_data_t;

typedef struct {  // yak parameters
	int32_t bf_shift, bf_n_hash;
	int32_t k, w, is_HPC;  // w==1 is "enumerate all kmers", 0 for "minimizers only"
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
} yak_copt_t;

typedef struct { // blocked bloom filter
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;

typedef struct {  // count hashtable
    // ?????? how does khashl_set's `key` value work though?
	yak_ct_t *h;  // it's a khashl_set hashtable; kmer with lower bits all zero is hashed; lower bits of the kh_key(h, hash_value) is used as the counter
	yak_bf_t *b;
} ha_ct1_t;

typedef struct {  // count hashtable
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	ha_ct1_t *h;  // ensembled hashtables
} ha_ct_t;

typedef struct {  // generate histogram
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {  // generate histogram
	const ha_ct_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

typedef struct {  // shrink a hashtable
	int min, max;
	ha_ct_t *h;
} shrink_aux_t;

typedef struct { // positional hash table
	yak_pt_t *h;  // an ensemble of hashtables
	uint64_t n;
	ha_idxpos_t *a;
} ha_pt1_t;

struct ha_pt_s { // positional hash table
	int k, pre;
	uint64_t tot, tot_pos;  // tot is set to be the tot in corresponding ct table; 
	ha_pt1_t *h;  // ensembled
};

typedef struct { // positional hash table
	const ha_ct_t *ct;
	ha_pt_t *pt;
} pt_gen_aux_t;

//************************
//       anchor
//************************

typedef struct { // this struct is not strictly necessary; we can use k_mer_pos instead, with modifications
	uint64_t srt;  // compressed for soring. upper bits is read ID + strand, like so: an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->other_off;
	uint32_t self_off:31, good:1;
	uint32_t other_off; // pos of the minimizer on the other read
} anchor1_t;

typedef struct {
	int n, good;
	const ha_idxpos_t *a;  // an array of minimize positional info (supplied by `ha_idx`)
} seed1_t;

struct ha_abuf_s {
	uint64_t n_a, m_a;  // n_a: total kmer occurrence counts
	uint32_t old_mz_m;
	ha_mz1_v mz;
	seed1_t *seed;
	anchor1_t *a;  // an array of anchor entries; will be radix-sorted
};

// #define an_key1(a) ((a).srt)
// #define an_key2(a) ((a).self_off)
// KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)  // 8 is (byte) "sizeof_key" => 64bits
// KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)  // 4 is (byte) "sizeof_key" => 32bits

// typedef struct
// {
// 	uint32_t readID:30, strand:1, good:1;
// 	uint32_t offset, self_offset;
// } k_mer_hit;

//************************
//       Hash_Table
//************************

typedef struct
{
    uint32_t offset;
    uint32_t readID:31, rev:1;
} k_mer_pos;

typedef struct
{
    k_mer_pos* list;
    uint64_t length;
    uint64_t size;
    uint8_t direction;
    uint64_t end_pos;
} k_mer_pos_list;

typedef struct
{
    int C_L[CIGAR_MAX_LENGTH];
    char C_C[CIGAR_MAX_LENGTH];
    int length;
} CIGAR;

typedef struct
{
  ///the begining and end of a window, instead of the whole overlap
  uint64_t x_start;
  uint64_t x_end;
  int y_end;
  int y_start;
  int extra_begin;
  int extra_end;
  int error_threshold;
  int error;
  CIGAR cigar;
} window_list;

typedef struct
{
    window_list* buffer;
    int32_t length;
    int32_t size;
} window_list_alloc;

typedef struct
{
    uint64_t* buffer;
	uint32_t length;
	uint32_t size;
} Fake_Cigar;

typedef struct
{
    uint32_t x_id;
    ///the begining and end of the whole overlap
    uint32_t x_pos_s;
    uint32_t x_pos_e;
    uint32_t x_pos_strand;

    uint32_t y_id;
    uint32_t y_pos_s;
    uint32_t y_pos_e;
    uint32_t y_pos_strand;

    uint32_t overlapLen;
    int32_t shared_seed;  // max_score; used by ks_introsort
    uint32_t align_length;
    uint8_t is_match;
    uint8_t without_large_indel;
    int8_t strong;
    uint32_t non_homopolymer_errors;

    window_list* w_list;
    uint32_t w_list_size;
    uint32_t w_list_length;
    Fake_Cigar f_cigar;

    window_list_alloc boundary_cigars;
} overlap_region;

typedef struct
{
    overlap_region* list;
    uint64_t size;
    uint64_t length;
    int64_t mapped_overlaps_length;
} overlap_region_alloc;

typedef struct
{
	uint32_t readID:30, strand:1, good:1;
	uint32_t offset, self_offset;
} k_mer_hit;

typedef struct {
	int32_t *score;
	int64_t *pre;
	int32_t *indels;
	int32_t *self_length;
	int64_t *tmp; // MUST BE 64-bit integer
	int64_t length;
	int64_t size;
} Chain_Data;

typedef struct
{
    k_mer_hit* list;  // a sorted list of minimizer hits
    long long length;  // total minizers (ab->n_a)
    long long size;  // capacity (ab->m_a)
    Chain_Data chainDP;
} Candidates_list;

//************************
//       Correct
//************************

typedef struct
{
    long long read_length;
    long long window_length;
    long long window_num;
    long long window_start;
    long long window_end;
    long long tail_length;
    int terminal;
}Window_Pool;

typedef struct
{
    /**[0-1] bits are type:**/
    /**[2-31] bits are length**/
    char current_operation;
    int current_operation_length;
    uint32_t* record;
    uint64_t size;
    uint64_t length;  // i think it's practically out degree of the corresponding vertex
    uint32_t new_read_length;


    char* lost_base;
    uint64_t lost_base_size;
    uint64_t lost_base_length;
    

}Cigar_record;


typedef struct
{
  long long length;
  long long size;
  Cigar_record* buffer;
}Cigar_record_alloc;


typedef struct
{
    ////the position of snp in read itself
    uint32_t site;
    ////the overlapID
    uint32_t overlapID;
    ////the position of snp in that overlap
    uint32_t overlapSite;
    ///there are several types: 0: equal to read 1: not equal to read, but it is a mismatch 2: is a gap
    uint8_t type;
    ///misbase
    char misBase;
}haplotype_evdience;



typedef struct
{
    ///the id of this snp
    uint32_t id;
    uint32_t overlap_num;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t occ_2;
    uint32_t homopolymer_num;
    uint32_t non_homopolymer_num;
    int score;
    ////the position of snp in read itself
    uint32_t site;
    uint8_t is_homopolymer;
}
SnpStats;



typedef struct
{
    uint32_t beg;
    uint32_t end;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t homopolymer_num;
    uint32_t non_homopolymer_num;
    uint32_t is_remove;
}
Snp_ID_Vector;


typedef struct
{
    long long IDs_size;
    long long IDs_length;
    long long max_snp_id;
    Snp_ID_Vector* IDs;

    long long buffer_size;
    long long buffer_length;
    uint32_t* buffer;
}
Snp_ID_Vector_Alloc;

#define Get_DP_Backtrack_Column(matrix, i) (matrix.backtrack + matrix.snp_num * i)
#define Get_DP_Backtrack_Column_Length(matrix, i) (matrix.snp_num)


typedef struct
{
    // uint32_t snp_size;
    // uint32_t snp_num;
    // uint32_t* max;
    // uint32_t* colum_len;
    // uint32_t* colum;
    // uint32_t matrix_size;

    uint32_t snp_num;
    
    uint8_t* visit;
    uint32_t* max;
    uint64_t* max_for_sort;



    uint32_t snp_size;
    
    uint32_t* backtrack_length;
    uint32_t* backtrack;
    uint32_t backtrack_size;


    uint32_t* buffer;
    uint32_t* max_buffer;


    int max_snp_num;
    ///int max_snp_ID;
    int max_score;
    int current_snp_num;


    Snp_ID_Vector_Alloc SNP_IDs;
}
DP_matrix;


#define Get_SNP_Martix_Size(matrix) (matrix.snp * matrix.overlap)
#define Get_SNP_Vector(matrix, i) (matrix.snp_matrix + matrix.overlap * i)
#define Get_SNP_Vector_Length(matrix) (matrix.overlap)
#define Get_Result_SNP_Vector(matrix) (matrix.snp_matrix + matrix.overlap*matrix.snp)

typedef struct
{
    haplotype_evdience* list;
    uint32_t sub_list_start;
    uint32_t sub_list_length;
    uint32_t length;
    uint32_t size;
    
    /****************************may have bugs********************************/
    uint8_t flag[WINDOW_MAX_SIZE];
    /****************************may have bugs********************************/

    uint32_t available_snp;
    uint32_t core_snp;
    uint32_t snp;
    uint32_t overlap;
    int8_t* snp_matrix;
    uint32_t snp_matrix_size;
    SnpStats* snp_stat; 
    SnpStats result_stat;
    uint32_t snp_stat_size;

    DP_matrix dp;
}
haplotype_evdience_alloc;

typedef struct
{
    char* corrected_read;
    long long corrected_read_length;
    long long last_boundary_length;
    long long corrected_read_size;
    long long corrected_base;

    uint64_t* overlapID;
    uint64_t length;
    uint64_t lengthNT;
    uint64_t size;
    uint64_t start_i;

    /****************************may have bugs********************************/
    // char overlap_region[WINDOW + THRESHOLD*2 + 10];
    // char overlap_region_group[GROUP_SIZE][WINDOW + THRESHOLD*2 + 10];
    // char path[WINDOW + THRESHOLD*2 + 10];
    // Word matrix_bit[((WINDOW + 10)<<3)];

    char overlap_region[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char overlap_region_group[GROUP_SIZE][WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char path[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];

    char path_fix[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    char overlap_region_fix[WINDOW_MAX_SIZE + THRESHOLD_MAX_SIZE*2 + 10];
    Word matrix_bit[((WINDOW_MAX_SIZE + 10)<<3)];
    /****************************may have bugs********************************/

    int path_length;
    __m128i Peq_SSE[256];
} Correct_dumy;

typedef struct
{
    Correct_dumy dumy;
    Cigar_record cigar;
    Cigar_record tmp_cigar;
    long long obtained_cigar_length;
}
Round2_alignment;

//************************
//       overlaps
//************************

typedef struct {
	int threadID;
    int thread_num;
    int check_cross;
    asg_t *g;
} para_for_simple_bub;

typedef struct { ///query is the read itself
	uint64_t qns;
	uint32_t qe, tn, ts, te;
	uint32_t ml:31, rev:1;
	uint32_t bl:31, del:1;  // bl is set to qe-qs?
	uint8_t el;
	uint8_t no_l_indel;
} ma_hit_t;  // equals to one paf line

typedef struct {
	ma_hit_t* buffer;
    uint32_t size;
    uint32_t length;
	uint8_t is_fully_corrected;
	uint8_t is_abnormal;
} ma_hit_t_alloc;  // paf

typedef struct {
	uint32_t s:31, del:1, e;
	uint8_t c;
} ma_sub_t;

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
	uint8_t strong;
	uint8_t el;
	uint8_t no_l_indel;
} asg_arc_t;


typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;



typedef struct {
	uint32_t len:31, del:1;
	uint8_t c;  // from coverage_cut
} asg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	asg_arc_t *arc;
	uint32_t m_seq, n_seq:31, is_symm:1;
	uint32_t r_seq;

	asg_seq_t *seq;
	uint64_t *idx;

	uint8_t* seq_vis;  // 1 if the node belongs to a bubble, 2 if it's a cross node (?)

	uint32_t n_F_seq;
	ma_utg_t* F_seq;
} asg_t;  // built by ma_sg_gen

typedef struct { size_t n, m; uint64_t *a; } asg64_v;


typedef struct { size_t n, m; ma_utg_t *a; } ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
} ma_ug_t;

typedef struct {
	uint32_t utg:31, ori:1, start, len;
} utg_intv_t;

/******************
 * Bubble popping *
 ******************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of positive reads
	uint32_t m; // max count of negative reads
	uint32_t np; // max count of non-positive reads
	uint32_t nc; // max count of reads, no matter positive or negative
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
	//s: state, s=0, this edge has not been visited, otherwise, s=1
} binfo_t;

typedef struct {
	///all information for each node
	binfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;


typedef struct {  // used for removing simple circle unitig
	kvec_t(uint64_t) Nodes; 
	kvec_t(uint64_t) Edges; 
	uint32_t pre_n_seq, seqID; 
} C_graph;

typedef struct {
	kvec_t(uint8_t) a;
	uint32_t i;
} kvec_t_u8_warp;

typedef struct {
	kvec_t(uint32_t) a;
	uint32_t i;
} kvec_t_u32_warp;

typedef struct {
	kvec_t(int32_t) a;
	uint32_t i;
} kvec_t_i32_warp;

typedef struct {
	kvec_t(uint64_t) a;
	uint64_t i;
} kvec_t_u64_warp;

typedef struct {
	kvec_t(asg_arc_t) a;
	uint64_t i;
}kvec_asg_arc_t_warp;

void sort_kvec_t_u64_warp(kvec_t_u64_warp* u_vecs, uint32_t is_descend);


typedef struct {
	uint32_t q_pos;
	uint32_t t_pos;
	uint32_t t_id;
	uint32_t is_color;
} Hap_Align;

typedef struct {
	kvec_t(Hap_Align) x;
	uint64_t i;
} Hap_Align_warp;

typedef struct {
	buf_t* b_0;
    uint32_t untigI;
    uint32_t readI;
    uint32_t offset;
} rIdContig;

typedef struct {
	uint64_t len;
	uint32_t* index;
} R_to_U;

typedef struct {
	asg_t* g;

	asg_arc_t *av;
	uint32_t nv;
	uint32_t av_i;

	asg_arc_t* new_edges;
	uint32_t new_edges_n;
	uint32_t new_edges_i;
} Edge_iter;

typedef struct {
	uint32_t father_occ;
	uint32_t mother_occ;
	uint32_t ambig_occ;
	uint32_t drop_occ;
	uint32_t total;
} Trio_counter;

//************************
//         POA
//************************
typedef struct
{
    uint64_t in_node;
    uint64_t out_node;
    ///0 is match，1 is mismatch，2 means y has more bases, 3 means x has more bases
    uint64_t weight;
    uint64_t num_insertions;
    uint64_t length;
    uint64_t self_edge_ID;
    uint64_t reverse_edge_ID;
} Edge;

typedef struct
{
    Edge* list;
    uint64_t size;
    uint64_t length;
    uint64_t delete_length;
} Edge_alloc;

typedef struct
{
    long long index;
} RSet;

typedef struct
{
    uint64_t ID;
    uint64_t weight;
    ///number of deletion end with current node
    uint64_t num_insertions;
    char base;
    Edge_alloc mismatch_edges;
    Edge_alloc deletion_edges;
    Edge_alloc insertion_edges;

} Node;

typedef struct
{
    uint64_t* list;
    uint8_t* visit;
    uint64_t size;
    uint64_t length;

    uint64_t* iterative_buffer;
    uint8_t* iterative_buffer_visit;
    uint64_t iterative_i;
} topo_Sorting_buffer;

typedef struct
{
    ///has a indivial start node 0
    Node* list;
    topo_Sorting_buffer sort;
    uint64_t size;
    uint64_t length;
    uint64_t delete_length;
} Node_alloc;

typedef struct
{
    uint64_t g_n_nodes;
    uint64_t g_n_edges;
    uint64_t g_next_nodeID;
    Node_alloc g_nodes;

    Queue node_q;
    char* seq;
    uint64_t seqID;
    uint64_t s_start_nodeID;
    uint64_t s_end_nodeID;
} Graph;