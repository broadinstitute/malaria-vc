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

    call tasks.IndexFasta {
        input:
            ref_fasta = ref_fasta
    }

    scatter (reads_bam in reads_bams) {
        call tasks.AlignSortDedupReads {
            input:
                reads_bam    = reads_bam,
                ref_fasta    = ref_fasta,
                ref_idx_dict = IndexFasta.ref_idx_dict,
                ref_idx_fai  = IndexFasta.ref_idx_fai,
                ref_idx_amb  = IndexFasta.ref_idx_amb,
                ref_idx_ann  = IndexFasta.ref_idx_ann,
                ref_idx_bwt  = IndexFasta.ref_idx_bwt,
                ref_idx_pac  = IndexFasta.ref_idx_pac,
                ref_idx_sa   = IndexFasta.ref_idx_sa
        }

        call tasks.BaseRecalibrator_1 {
            input:
                aligned_bam     = AlignSortDedupReads.aligned_bam,
                aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
                ref_fasta       = ref_fasta, 
                ref_idx_dict    = IndexFasta.ref_idx_dict, 
                ref_idx_fai     = IndexFasta.ref_idx_fai
        }

        call tasks.BaseRecalibrator_2 {
            input:
                bqsr_table      = BaseRecalibrator_1.table,
                aligned_bam     = AlignSortDedupReads.aligned_bam,
                aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx,
                ref_fasta       = ref_fasta, 
                ref_idx_dict    = IndexFasta.ref_idx_dict, 
                ref_idx_fai     = IndexFasta.ref_idx_fai
        }

        call tasks.HaplotypeCaller {
            input:
                aligned_bam     = AlignSortDedupReads.aligned_bam,
                aligned_bam_idx = AlignSortDedupReads.aligned_bam_idx, 
                bqsr_table      = BaseRecalibrator_1.table,
                ref_fasta       = ref_fasta, 
                ref_idx_dict    = IndexFasta.ref_idx_dict,
                ref_idx_fai     = IndexFasta.ref_idx_fai
        }
    }

    call GenotypeGVCFs {
        input:
            vcf_files    = HaplotypeCaller.vcf,
            ref_fasta    = ref_fasta,
            ref_idx_dict = IndexFasta.ref_idx_dict,
            ref_idx_fai  = IndexFasta.ref_idx_fai
    }
    
    call HardFiltration { 
        input: 
            vcf          = GenotypeGVCFs.out_vcf,
            vcf_tbi      = GenotypeGVCFs.out_vcf_tbi,
            ref_fasta    = ref_fasta,
            ref_idx_dict = IndexFasta.ref_idx_dict,
            ref_idx_fai  = IndexFasta.ref_idx_fai
    }

    call SnpEff {
        input:
            vcf       = HardFiltration.out_vcf,
            vcf_tbi   = HardFiltration.out_vcf_tbi,
            ref_fasta = ref_fasta
    }

}

