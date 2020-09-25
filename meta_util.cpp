#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "meta_util.h"
#include <assert.h>
#include "Process_Read.h"
#include "kthread.h"
#include "CommandLines.h"  // for Get_T()
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>  // debug, for printing uint64
#define MIN(a,b) ((a)<=(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

///////////////* moving average * ///////////////
// improve: merge into the relevant function and use queue
void moving_averagel(const uint16_t *counts, uint32_t l, uint16_t w, uint16_t *buf){
    uint16_t i0=0;
    uint16_t i1=w-1;
    double tmp = 0;
    for (uint16_t i=0; i<w; i++) {
        tmp+= (double)counts[i]/w;  // not the best way
    }
    uint16_t i;
    for (i=0; i<l-w+1; i++){
        buf[i] = (uint16_t)tmp;
        i1++;
        tmp = tmp - (double)(counts[i0])/w + (double)(counts[i1])/w;
        i0++;
    }
}

double meanl(const uint16_t *counts, uint32_t l){
    double ret = 0;
    for (uint32_t i=0; i<l; i++){
        ret+=counts[i];
    }
    ret = ret/l;
    if (ret==0) ret = 1;  // debug
    // overflow-proof version:
    // uint32_t t = 1;
    // for (uint32_t i=0; i<l; i++){
    //     ret+=((double)counts[i]-ret)/t;
    //     t++;
    // }
    return ret;
}

uint64_t medianl(const uint16_t *counts, uint32_t l){
    return counts[l/2];
}

double stdl(const uint16_t *counts, uint32_t l, double mean){
    // note: should not overflow...
    if (l==0) return (double)-1;
    double ret = 0;
    for (uint32_t i=0; i<l; i++){
        if (counts[i]==0) ret+=1; // debug
        ret+=pow(counts[i]-mean, 2);
    }
    if (ret==0) return (double)0;
    return sqrt(ret/l);
}

// uint32_t saturationl(const uint64_t *counts, uint32_t l, uint64_t threshold){
//     uint16_t ret = 0;
//     if (uint32_t i=0; i<l; i++){
//         if (counts[i]>threshold) ret++;
//     }
//     return ret;
// }

///////////////* probing *//////////////////
// hamt_peaks_t hamt_find_peaks(uint32_t n_cnt, const uint64_t *cnt, uint16_t l_smooth){
//     hamt_peaks_t peaks = {0, 0, 0xffffffffffffffff, 0};
//     // collect smoothed peak values
//     uint16_t *buf = (uint16_t*) malloc(sizeof(uint16_t)*(n_cnt-l_smooth+1));
//     moving_averagel(cnt, n_cnt, l_smooth, buf);
//     double acc = 0;
//     for (uint32_t i=1; i<n_cnt-l_smooth; i++){
//         if (buf[i]>buf[i-1] && buf[i]>buf[i+1]){
//             printf(" - peak at %d, ", i);
//             printf("%" PRIu64 "\n", buf[i]);
//             peaks.n++;
//             if (buf[i]>peaks.max) peaks.max = buf[i];
//             if (buf[i]<peaks.min) peaks.min = buf[i];
//             acc += (double)buf[i]/500;
//         }
//     }
//     peaks.mean = acc/peaks.n*500;
//     free(buf);
//     return peaks;
// }

//////////////////* marking *////////////////////
#define HAMT_READCOV_LOW 0x1
#define HAMT_READCOV_REASONABLE 0x2
#define HAMT_READCOV_HIGH 0x4
#define HAMT_READDIV_UNIFORMED 0x10
#define HAMT_READDIV_LONGLOW 0x20
#define HAMT_READDIV_REASONABLE 0x40
#define HAMT_READDIV_HIGH 0x80

#define HAMT_DISCARD 0x1
#define HAMT_VIA_MEDIAN 0x2
#define HAMT_VIA_LONGLOW 0x4
#define HAMT_VIA_KMER 0x8

uint8_t decide_category(double mean, double std, uint16_t *buf, uint32_t l){  // used for the initial pass
    // buf is UNSORTED kmer count along the read, and l is the length of buf
    uint8_t flag = 0;
    int low_threshold = asm_opt.diginorm_coverage/2 >50? 50 : asm_opt.diginorm_coverage/2;
    // coverage
    if (mean<low_threshold) flag|=HAMT_READCOV_LOW;
    else if (mean<asm_opt.diginorm_coverage) flag|=HAMT_READCOV_REASONABLE;
    else flag|=HAMT_READCOV_HIGH;
    // divergence
    uint32_t cnt = 0;
    uint32_t max_cnt = 0;
    double flt_low = MAX(mean*0.1, 1); 
    double flt_break = MAX(mean*0.3, 1);
    if (mean>=low_threshold){  // read kmers are reasonable prevalent
        for (uint32_t i=0; i<l-1; i++){
            if (buf[i]<=flt_low)
                cnt+=1;
            else if (buf[i]>flt_break){
                if (cnt>max_cnt) max_cnt = cnt;
                cnt = 0;
            }
        }
        if ((double)max_cnt/l>0.3) {flag|=HAMT_READDIV_LONGLOW;}
    }
    if ((std/mean)<0.08 || std<20) flag|=HAMT_READDIV_UNIFORMED;  // very uniformed read (or just rare), 8% is roughly 8/255
    else if ((std/mean)<0.20) flag|=HAMT_READDIV_REASONABLE;  // regular, 20% is roughly 50/255
    else flag|=HAMT_READDIV_HIGH;
    return flag;
}

uint8_t decide_drop( double mean, double std, uint16_t runtime_median, int round, uint8_t initmark){
    uint8_t code = 0;
    if (runtime_median<asm_opt.diginorm_coverage) code = HAMT_VIA_MEDIAN;  // note to self: 50 or 100, this only reduce like less than 50% reads in real data.
    else if (runtime_median<asm_opt.diginorm_coverage*2 and (initmark&HAMT_READDIV_LONGLOW)) code= HAMT_VIA_LONGLOW;  // keep long low-coverage read
    else if (1) code = HAMT_DISCARD;
    return code;
}

//////////////////* del overlaps *////////////////////
typedef struct{
    All_reads *rs;
    int which_paf;
    uint64_t *cnt;

    int diff_abs;
    double diff_fold;
} ktforbuf_del_by_cov_t;

static void worker_hamt_del_ovlp_by_coverage(void* data, long idx_for, int tid){
    ktforbuf_del_by_cov_t *s = (ktforbuf_del_by_cov_t*) data;
    ma_hit_t *h = NULL;
    ma_hit_t_alloc *x = s->which_paf==0? &(s->rs->paf[idx_for]) : &(s->rs->reverse_paf[idx_for]);
    uint8_t code_query = s->rs->mask_readtype[idx_for];  // diginorm_coverage bit flag of the query read
    int median_query = s->rs->median[idx_for];
    median_query = median_query!=0? median_query : 1;  // just in case for divided by 0
    uint8_t code_target = 0;
    int median_target = 0;

    int diff_abs;
    double diff_fold;
    int tmp_cnt = 0;

    for (int i = 0; i < x->length; i++){
        h = &(x->buffer[i]);  // the overlap
        if(h->del) continue;
        code_target = s->rs->mask_readtype[h->tn];
        assert(code_target!=0);
        if ((code_query & HAMT_READCOV_LOW && code_target & HAMT_READCOV_HIGH) || (code_target & HAMT_READCOV_LOW && code_query & HAMT_READCOV_HIGH)){
            h->del = 1;
            tmp_cnt+=1;
        }else{
            median_target = s->rs->median[h->tn];
            median_target = median_target!=0? median_target : 1;  // just in case for divided by 0
            if (median_query>median_target){
                diff_abs = median_query - median_target;
                diff_fold = (double)median_query/median_target;
            }else{
                diff_abs = median_target - median_query;
                diff_fold = (double)median_target/median_query;
            }
            if (diff_abs >= s->diff_abs) {h->del = 1; tmp_cnt+=1; continue;}
            if ((median_query>=8 || median_target >=8) && (diff_fold >= s->diff_fold)) {h->del = 1; tmp_cnt+=1; continue;}
        }
    }

    *s->cnt += tmp_cnt;   // not safe but ok?
}

void hamt_del_ovlp_by_coverage(All_reads *rs, hifiasm_opt_t asm_opt, int diff_abs, double diff_fold){
    double startTime = Get_T();

    uint64_t cnt_dropped_ovlp_p = 0;
    uint64_t cnt_dropped_ovlp_r = 0;
    ktforbuf_del_by_cov_t s;
    s.rs = rs;
    s.diff_abs = diff_abs;
    s.diff_fold = diff_fold;

    s.which_paf = 0;  // use rs->paf
    s.cnt = &cnt_dropped_ovlp_p;  // sancheck counter
    kt_for(asm_opt.thread_num, worker_hamt_del_ovlp_by_coverage, &s, rs->total_reads);
    
    if (!asm_opt.is_disable_phasing){
        s.which_paf = 1;  // use rs->reverse_paf
        s.cnt = &cnt_dropped_ovlp_r;
        kt_for(asm_opt.thread_num, worker_hamt_del_ovlp_by_coverage, &s, rs->total_reads);
    }
    fprintf(stderr, "[meta_util::%s] took %0.2fs, masked %" PRIu64 " overlaps for paf, %" PRIu64 " for reverse_paf.\n", __func__, Get_T()-startTime, cnt_dropped_ovlp_p, cnt_dropped_ovlp_r);

}