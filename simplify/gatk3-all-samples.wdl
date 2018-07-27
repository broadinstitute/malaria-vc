# 
#  gatk3-all-samples.wdl: 
#     - standardized pipeline for calling snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

import "tasks.wdl" as tasks

workflow GATK3_Joint_Call_and_Annotate {

    Array[File]+ reads_bams
    File         ref_fasta

    Array[File]? bqsr_known_sites  # empty array will skip BQSR

    scatter (reads_bam in reads_bams) {
        call tasks.AlignSortDedupReads {
            input:
                reads_bam = reads_bam,
                ref_fasta = ref_fasta
        }

        call tasks.BaseRecalibrator {
            input:
                aligned_bam     = AlignSortDedupReads.aligned_bam,
                aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
                known_sites     = bqsr_known_sites,
                ref_fasta       = ref_fasta
        }

        call tasks.HaplotypeCaller {
            input:
                aligned_bam     = BaseRecalibrator.out_bam,
                aligned_bam_idx = BaseRecalibrator.out_bam_idx, 
                ref_fasta       = ref_fasta
        }
    }

    call tasks.GenotypeGVCFs {
        input:
            vcf_files    = HaplotypeCaller.vcf,
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

