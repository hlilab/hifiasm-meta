#include <stdint.h>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>  // debug, for printing uint64
#include <time.h>

#include "htab.h"
#include "CommandLines.h"


// copynpaste from htab.cpp
#define HAF_COUNT_EXACT  0x1
#define HAF_COUNT_ALL    0x2
#define HAF_RS_WRITE_LEN 0x4
#define HAF_RS_WRITE_SEQ 0x8
#define HAF_RS_READ      0x10
#define HAF_CREATE_NEW   0x20
#define HAMTF_FORCE_DONT_INIT 0x40  // meta: override HAF_RS_WRITE_LEN to not allow init_all_reads of ha_count


int hamt_rkp(void){
	hamt_read_kmer_profile();
}