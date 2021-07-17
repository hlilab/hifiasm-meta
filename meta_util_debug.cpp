#include <stdio.h>
#include "Process_Read.h"
#include "CommandLines.h"
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>

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
void hamt_dump_selected_read_names(hifiasm_opt_t *asm_opt, All_reads *rs){
    char *base_name = (char*)malloc(strlen(asm_opt->bin_base_name)+25);
    sprintf(base_name, "%s.readlist", asm_opt->bin_base_name);
    FILE *fp = fopen(base_name, "w");
    int tot = 0;
    for (uint64_t i=0; i<rs->total_reads; i++){
        if (rs->mask_readnorm[i]&1){continue;}
        fprintf(fp, "%.*s\n", (int)Get_NAME_LENGTH((*rs), i), Get_NAME((*rs), i));
        tot++;
    }
    fprintf(stderr, "[M::%s] wrote %d / %d read names\n", __func__, tot, (int)rs->total_reads);
    fclose(fp);
    free(base_name);

}


void hamt_dump_ovec_read_error_count_and_kmerinfo(hifiasm_opt_t *asm_opt, All_reads *rs){
    // (only for bin files since r7-ish)
    // output the total number of bases corrected for each read
    char *bin_base_name = (char*)malloc(strlen(asm_opt->bin_base_name)+25);
    sprintf(bin_base_name, "%s.ec", asm_opt->bin_base_name);
    
    char *output_file_name = (char*)malloc(strlen(bin_base_name)+25);
    
    if (strcmp(asm_opt->output_file_name, "hifiasm.asm")==0){
        fprintf(stderr, "[W::%s] -o not given, defaulting to -B's value (%s).\n", __func__, asm_opt->bin_base_name);
        sprintf(output_file_name, "%s.ovec_count_dump", bin_base_name);
    }else{
        sprintf(output_file_name, "%s.ovec_count_dump", asm_opt->output_file_name);
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
        fprintf(fp, "%s\t%d\t%f\t%d\n", readname_s, (int)rs->median[i], rs->std[i], (int)rs->nb_error_corrected[i]);
    }

    fclose(fp);
    free(output_file_name);
    free(readname_s);
    destory_All_reads(rs);

}

void hist_readlength(All_reads *rs){
    int bins, max_length=50000, min_length = 1000, step=500;
    bins = max_length/step;
    int under=0, over=0;
    int topping = 50;
    int d = 5000, topped;  // discrete bar to ascii
    
    int l, x, last_idx=0;
    uint32_t buf[bins];
    memset(buf, 0, bins*sizeof(uint32_t));

    for (uint64_t i=0; i<rs->total_reads; i++){
        l = (int) (Get_READ_LENGTH((*rs), i));
        if (l<min_length){
            under++;
            continue;
        }else if (l>=max_length){
            over++;
            continue;
        }
        l = (l-min_length)/step;
        buf[l]++;
    }
    for (int i=bins-1; i>=0; i--){
        if (buf[i]>0){
            last_idx = i;
            break;
        }
    }

    // hist
    if (under==0){
        fprintf(stderr, "[M::%s] <%5.1fk: 0\n", __func__, (float)min_length/1000);
    }else{
        fprintf(stderr, "[M::%s] <%5.1fk: ", __func__, (float)min_length/1000);
        topped = 0;
        x = under/d;
        if (x>topping) {topped = 1;}
        else{topped = 0;}
        for (int j=0; j<x; j++){
            fputc(']', stderr);
            if (j>topping) break;
        }
        fprintf(stderr, " %d", x);
        if (topped){fputc(')', stderr);}
        fputc('\n', stderr);
    }
    for (int i=0; i<last_idx+1; i++){
        if (buf[i]==0){
            fprintf(stderr, "[M::%s] %.1fk: 0\n", __func__, (float)(min_length+step*i)/1000);
        }else{
            fprintf(stderr, "[M::%s] %.1fk: ", __func__, (float)(min_length+step*i)/1000);
            if (buf[i]==0){  // this bin doesn't have any read
                fputc('0', stderr);
            }else{
                x = buf[i]/d;
                if (x>topping) {topped = 1;}
                else{topped = 0;}
                for (int j=0; j<x+1; j++){
                    fputc(']', stderr);
                    if (j>topping) break;
                }
                if (topped){fputc(')', stderr);}
                fprintf(stderr, " %" PRIu32 "", buf[i]);
            }
            fputc('\n', stderr);
        }
    }
    if (over==0){
        fprintf(stderr, "[M::%s] >%.1fk: 0\n", __func__, ((float)max_length)/1000);
    }else{
        fprintf(stderr, "[M::%s] >%.1fk: ", __func__, ((float)max_length)/1000);
        x = over/d;
        if (x>100) {topped = 1;}
        else{topped = 0;}
        for (int j=0; j<x; j++){
            fputc(']', stderr);
        }
        if (topped){fputc('+', stderr);}
        fputc('\n', stderr);
    }

}