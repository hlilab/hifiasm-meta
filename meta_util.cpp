#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "meta_util.h"
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

uint8_t decide_category(double mean, double std, uint16_t *buf, uint32_t l){  // used for the initial pass
    // buf is UNSORTED kmer count along the read, and l is the length of buf
    uint8_t flag = 0;
    // coverage
    if (mean<30) flag|=HAMT_READCOV_LOW;
    else if (mean<100) flag|=HAMT_READCOV_REASONABLE;
    else flag|=HAMT_READCOV_HIGH;
    // divergence
    uint32_t cnt = 0;
    uint32_t max_cnt = 0;
    double flt_low = MAX(mean*0.1, 1); 
    double flt_break = MAX(mean*0.3, 1);
    if (mean>=30){  // read kmers are reasonable prevalent
        for (uint32_t i=0; i<l-1; i++){
            // buf_diff[i] = (double)(buf[i+1]-buf[i])/buf[i+1];
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

int decide_drop( double mean, double std, uint16_t runtime_median, int round, uint8_t initmark){
    if (initmark==0) return 0;
	if (round==1){
        if (runtime_median<50) return 0;  // note to self: 50 or 100, this only reduce like less than 50% reads in real data.
        if (runtime_median<50 and (initmark&HAMT_READDIV_LONGLOW)) return 0;  // keep long low-coverage read
    }
    return 1;
}

