# 
#  GATK3-merge-samples.wdl: 
#     - standardized pipeline for merging per-sample HaplotypeCaller outputs,
#       calling SNPs and indels jointly, hard filtering, and annotating. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow GATK3_Multi_Sample {

    File ref_fasta
    File ref_idx_dict
    File ref_idx_fai

    call GenotypeGVCFs {
        input:
            ref_fasta    = ref_fasta,
            ref_idx_dict = ref_idx_dict,
            ref_idx_fai  = ref_idx_fai
    }
    
    call HardFiltration { 
        input: 
            vcf          = GenotypeGVCFs.out_vcf,
            vcf_tbi      = GenotypeGVCFs.out_vcf_tbi,
            ref_fasta    = ref_fasta,
            ref_idx_dict = ref_idx_dict,
            ref_idx_fai  = ref_idx_fai
    }

    call SnpEff {
        input:
            vcf       = HardFiltration.out_vcf,
            vcf_tbi   = HardFiltration.out_vcf_tbi,
            ref_fasta = ref_fasta
    }
}

