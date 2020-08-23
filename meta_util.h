#ifndef __META_UTILITY__
#define __META_UTILITY__
#include <stdint.h>

typedef struct hamt_peaks_s{
    uint64_t n, max, min, mean;
} hamt_peaks_t;

double meanl(const uint16_t *counts, uint32_t l);
double stdl(const uint16_t *counts, uint32_t l, double mean);
uint64_t medianl(const uint16_t *counts, uint32_t l);
void moving_averagel(uint16_t *counts, uint32_t l, uint16_t w, uint64_t *buf);
// hamt_peaks_t hamt_find_peaks(uint32_t n_cnt, const int64_t *cnt, uint16_t l_smooth);
uint8_t decide_category(double mean, double std, uint16_t *buf, uint32_t l);
int decide_drop(double mean, double std, uint16_t runtime_median, int round, uint8_t initmark);

#endif  //__META_UTILITY__