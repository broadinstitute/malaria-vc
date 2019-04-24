# 
#  gatk3-one-sample.wdl: 
#     - standardized pipeline for aligning reads and calling SNPs and
#       indels for a single sample using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow Align_HaplotypeCaller_One_Sample {

    File ref_fasta

    call tasks.AlignSortDedupReads {
        input:
            ref_fasta    = ref_fasta
    }

    call tasks.BaseRecalibrator {
        input:
            aligned_bam     = AlignSortDedupReads.aligned_bam,
            aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
            ref_fasta       = ref_fasta
    }

    call tasks.HaplotypeCaller {
        input:
            aligned_bam     = BaseRecalibrator.out_bam,
            aligned_bam_idx = BaseRecalibrator.out_bam_idx, 
            ref_fasta       = ref_fasta
    }
}
