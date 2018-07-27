# 
#  gatk3-merge-samples.wdl: 
#     - standardized pipeline for merging per-sample HaplotypeCaller outputs,
#       calling SNPs and indels jointly, hard filtering, and annotating. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow Merge_gVCFs_Genotype_Annotate {

    File ref_fasta

    call tasks.GenotypeGVCFs {
        input:
            ref_fasta    = ref_fasta
    }
    
    call tasks.HardFiltration { 
        input: 
            vcf          = GenotypeGVCFs.out_vcf,
            vcf_tbi      = GenotypeGVCFs.out_vcf_tbi,
            ref_fasta    = ref_fasta
    }

    call tasks.SnpEff {
        input:
            vcf       = HardFiltration.out_vcf,
            vcf_tbi   = HardFiltration.out_vcf_tbi,
            ref_fasta = ref_fasta
    }
}

