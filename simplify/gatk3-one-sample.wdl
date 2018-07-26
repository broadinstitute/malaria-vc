# 
#  GATK3_One_Sample.wdl: 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow GATK3_One_Sample {
    File ref_fasta
    File ref_idx_dict
    File ref_idx_fai

    call tasks.AlignSortDedupReads {
        input:
            ref_fasta = ref_fasta,
            ref_idx_dict = ref_idx_dict,
            ref_idx_fai = ref_idx_fai
    }

    call tasks.BQSR {
        input:
            aligned_bam = AlignSortDedupReads.aligned_bam,
            aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
            ref_fasta = ref_fasta, 
            ref_idx_dict = ref_idx_dict, 
            ref_idx_fai = ref_idx_fai
    }

    call tasks.HaplotypeCaller {
        input:
            aligned_bam = BQSR.out_bam,
            aligned_bam_idx = BQSR.out_bam_idx, 
            bqsr_table = BQSR.table,
            ref_fasta = ref_fasta, 
            ref_idx_dict = ref_idx_dict,
            ref_idx_fai = ref_idx_fai
    }
}
