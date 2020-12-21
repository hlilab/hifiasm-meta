#ifndef __META_U_DBG__
#define __META_U_DBG__
#include "CommandLines.h"
#include "Process_Read.h"

void hamt_dump_read_selection_mask_runtime(hifiasm_opt_t *asm_opt, All_reads *rs);
void hamt_dump_read_selection_mask(hifiasm_opt_t *asm_opt, All_reads *rs);
void hamt_dump_selected_read_names(hifiasm_opt_t *asm_opt, All_reads *rs);
void hamt_dump_ovec_read_error_count_and_kmerinfo(hifiasm_opt_t *asm_opt, All_reads *rs);
void hist_readlength(All_reads *rs);

#endif  // __META_U_DBG__