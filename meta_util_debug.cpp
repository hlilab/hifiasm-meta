#include <stdio.h>
#include "Process_Read.h"
#include "CommandLines.h"

void hamt_dump_read_selection_mask_runtime(hifiasm_opt_t *asm_opt, All_reads *rs){
    char *output_filename = (char*)malloc(strlen(asm_opt->output_file_name)+25);
    sprintf(output_filename, "%s.dump_read_mask", asm_opt->output_file_name);
    FILE *fp = fopen(output_filename, "w");

    char *readname_s = (char*)malloc(100);
    int readname_l;
    for (uint64_t i=0; i<rs->total_reads; i++){
        readname_l = rs->name_index[i+1]-rs->name_index[i];
        memcpy(readname_s, &rs->name[rs->name_index[i]], readname_l);
        readname_s[readname_l] = '\0';
        fprintf(fp, "%s\t%d\t%d\t%f\t%d\n", readname_s, (int)rs->mask_readnorm[i], (int)rs->median[i], rs->std[i], (int)rs->lowq[i]);
    }

    fclose(fp);
    free(output_filename);
    free(readname_s);
}


void hamt_dump_read_selection_mask(hifiasm_opt_t *asm_opt, All_reads *rs){
    // given bin files, write down the read selection mask
    char *bin_base_name = (char*)malloc(strlen(asm_opt->bin_base_name)+5);
    sprintf(bin_base_name, "%s.ec", asm_opt->bin_base_name);
    if (strcmp(bin_base_name, "hifiasm.asm")==0){
        fprintf(stderr, "[E::%s] must specify bin file via -B. Don't use -o, it's reserved for output name.\n", __func__);
        exit(1);
    }
    char *output_file_name = (char*)malloc(strlen(bin_base_name)+25);
    
    if (strcmp(asm_opt->output_file_name, "hifiasm.asm")==0){
        fprintf(stderr, "[W::%s] -o not given, defaulting to -B's value (%s).\n", __func__, asm_opt->bin_base_name);
        sprintf(output_file_name, "%s.readinfodump", bin_base_name);
    }else{
        sprintf(output_file_name, "%s.readinfodump", asm_opt->output_file_name);
    }

    fprintf(stderr, "[M::%s] loade from base name %s\n", __func__, bin_base_name);
    int ret = load_All_reads(rs, bin_base_name);
    if (!ret){
        fprintf(stderr, "[E::%s] can't load from bin files.\n", __func__);
        exit(1);
    }
    fprintf(stderr, "[M::%s] loaded all reads, dumping..\n", __func__);

    char *readname_s = (char*)malloc(100);
    int readname_l;
    FILE *fp = fopen(output_file_name, "w");
    for (uint64_t i=0; i<rs->total_reads; i++){
        readname_l = rs->name_index[i+1]-rs->name_index[i];
        memcpy(readname_s, &rs->name[rs->name_index[i]], readname_l);
        readname_s[readname_l] = '\0';
        fprintf(fp, "%s\t%d\t%d\t%f\n", readname_s, (int)rs->mask_readnorm[i], (int)rs->median[i], rs->std[i]);
    }

    fclose(fp);
    free(output_file_name);
    free(readname_s);
    destory_All_reads(rs);

}