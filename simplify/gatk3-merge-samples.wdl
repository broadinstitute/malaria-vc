# 
#  Malaria_WGS_GATK3.wdl: 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

## WORKFLOW DEFINITION
workflow GATK3_Merge_Samples {

    File ref_fasta
    File ref_idx_dict
    File ref_idx_fai

    # for base quality score recalibration
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] bqsr_intervals_to_exclude

    # variant quality control param
    # either "vqsr" or "hard_filtering"
    String variant_qc="hard_filtering"
    
    # vqsr params       
    # if variant_qc == "vqsr"
    # all of these params are required
    Float ts_filter_snp=99.0
    Float ts_filter_indel=99.5
    Int snp_max_gaussians=8
    Int indel_max_gaussians=4
    Int vqsr_mapping_qual_cap=70
    Array[String] snp_resources=["7g8_gb4,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/3d7_hb3.combined.final.vcf.gz"]
    Array[String] indel_resources=["7g8_gb4,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/3d7_hb3.combined.final.vcf.gz"]
    Array[String] snp_annotations=["QD", "FS", "SOR", "DP"]
    Array[String] indel_annotations=["QD", "FS"]

    call GenotypeGVCFs {
        input:
            ref_fasta = ref_fasta,
            ref_idx_dict = ref_idx_dict,
            ref_idx_fai = ref_idx_fai
    }
    
    # variant quality control
    if (variant_qc == "vqsr") { 
        # variant quality score recalibration
        # snp vqsr
        call VQSR as SnpVQSR {
            input:
                ref_fasta = ref_fasta,
                ref_idx_dict = ref_idx_dict,
                ref_idx_fai = ref_idx_fai
                intervals = interval_files_list,
                gvcf = GenotypeGVCFs.out,
                output_filename = "${run_name}.snp_vqsr.g.vcf",
                
                mode = "SNP",
                resources = snp_resources,
                resource_files = known_sites, 
                resource_file_indices = known_sites_indices,
                annotations = snp_annotations,
                ts_filter = ts_filter_snp,
                max_gaussians = snp_max_gaussians,
                mapping_qual_cap = vqsr_mapping_qual_cap
        }
        # indel vqsr
        call VQSR as IndelVQSR {
            input:
                ref_fasta = ref_fasta,
                ref_idx_dict = ref_idx_dict,
                ref_idx_fai = ref_idx_fai
                intervals = interval_files_list,
                gvcf = SnpVQSR.out,
                output_filename = "${run_name}.indel_vqsr.g.vcf",

                mode = "INDEL",
                resources = indel_resources,
                resource_files = known_sites, 
                resource_file_indices = known_sites_indices,
                annotations = indel_annotations,
                ts_filter = ts_filter_indel,
                max_gaussians = indel_max_gaussians,
                mapping_qual_cap = vqsr_mapping_qual_cap
        }
    } 
    if (variant_qc == "hard_filtering") { 
        call HardFiltration { 
            input: 
                vcf             = GenotypeGVCFs.out,
                ref_fasta       = ref_fasta,
                ref_idx_dict    = ref_idx_dict,
                ref_idx_fai     = ref_idx_fai
        }
    }

    call SnpEff {
        input:
            vcf = select_first([IndelVQSR.out, HardFiltration.out]),
            ref_fasta = ref_fasta, 
            output_filename = "${run_name}.snpeff.g.vcf"
    }
}

