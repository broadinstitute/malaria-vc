# 
#  gatk3-one-sample.wdl: 
#     - standardized pipeline for aligning reads and calling SNPs and
#       indels for a single sample using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow GATK3_One_Sample_HC {
    File ref_fasta
    File ref_idx_dict
    File ref_idx_fai

    call tasks.AlignSortDedupReads {
        input:
            ref_fasta    = ref_fasta
    }

    call tasks.BaseRecalibrator_1 {
        input:
            aligned_bam     = AlignSortDedupReads.aligned_bam,
            aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
            ref_fasta       = ref_fasta, 
            ref_idx_dict    = ref_idx_dict, 
            ref_idx_fai     = ref_idx_fai
    }

    call tasks.BaseRecalibrator_2 {
        input:
            bqsr_table      = BaseRecalibrator_1.table,
            aligned_bam     = AlignSortDedupReads.aligned_bam,
            aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
            ref_fasta       = ref_fasta, 
            ref_idx_dict    = ref_idx_dict, 
            ref_idx_fai     = ref_idx_fai
    }

    call tasks.HaplotypeCaller {
        input:
            aligned_bam     = AlignSortDedupReads.aligned_bam,
            aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx, 
            bqsr_table      = BaseRecalibrator_1.table,
            ref_fasta       = ref_fasta, 
            ref_idx_dict    = ref_idx_dict,
            ref_idx_fai     = ref_idx_fai
    }
}
