# 
#  Malaria_WGS_GATK3.wdl: 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

## WORKFLOW DEFINITION
workflow Malaria_WGS_GATK3 {

    ## config params
    # input data
    File ref			       # path to reference file
    File ref_dict
    File ref_index
    String run_name
    File sample_paths_file     # .tsv of (sample_name, sample_sam_path)
    File interval_files_list		   
    
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

    # hard filtering params
    # if variant_qc == "hard_filtering"
    # both of these params are required
    String snp_filter_expr="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    String indel_filter_expr="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
    
    # snpeff 
    File Pf3D7_gff
    Boolean use_snpeff=true


    ## task calls 
    # run pipeline on each sample, in parallel
    scatter(sample in read_tsv(sample_paths_file)) {
        String sample_name = sample[0]
        String sample_bam_path = sample[1]	

        call AlignSortDedupReads {
            input:
            sample_name = sample_name,
            reads_bam = sample_bam_path
        }

        # base quality score recalibration 
        call BQSR {
            input:
            ref_fasta = ref, 
            ref_dict = ref_dict, 
            ref_index = ref_index,
            bam = AlignSortDedupReads.aligned_bam,
            bam_index = AlignSortDedupReads.aligned_bam_index,
            sample_name = sample_name, 
            known_sites = known_sites,
            known_sites_indices = known_sites_indices,
            intervals_to_exclude = bqsr_intervals_to_exclude,
            output_table_name = sample_name + ".bqsr.table",
            output_bam_name = sample_name + ".bsqr.bam"
        }

        call HaplotypeCaller {
            input:
            ref = ref, 
            bam = BQSR.out,
            ref_dict = ref_dict,
            ref_index = ref_index, 
            bam_index = BQSR.out_index, 
            bqsr_table = BQSR.table,
            sample_name = sample_name, 
            intervals = interval_files_list
        }
    } # end scatter

    call GenotypeGVCFs {
        input:
        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index, 
        intervals = interval_files_list,
        vcf_files = HaplotypeCaller.vcf
    }
    
    # variant quality control
    if (variant_qc == "vqsr") { 
        # variant quality score recalibration
        # snp vqsr
        call VQSR as SnpVQSR {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index, 
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
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index, 
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
            ref = ref, 
            ref_dict = ref_dict, 
            ref_index = ref_index, 
            vcf = GenotypeGVCFs.out,
            snp_filter_expr = snp_filter_expr,
            indel_filter_expr = indel_filter_expr,
            output_filename = "${run_name}.hard_filtered.g.vcf"
        }
    }

    # add variant annotations using SnpEff
    File vcf = select_first([IndelVQSR.out, HardFiltration.out])
    if (use_snpeff == true) {
        call SnpEff {
            input:
            vcf = vcf,
            ref = ref, 
            Pf3D7_gff = Pf3D7_gff,
            output_filename = "${run_name}.snpeff.g.vcf"
        }
    }
    
    output { 
        File gvcf = select_first([SnpEff.out, IndelVQSR.out, HardFiltration.out])
    }
}

## TASK DEFINITIONS 

task AlignSortDedupReads {
    File reads_bam
    File ref_fasta
    File ref_idx_dict
    File ref_idx_fai
    File ref_idx_amb
    File ref_idx_ann
    File ref_idx_bwt
    File ref_idx_pac
    File ref_idx_sa

    command {
        set -ex -o pipefail
        _MEM="12G"
        _MEM_HALF="6G"

        # TODO: pull {read_group} from in_bam header
        bwa index ${ref_fasta}

        # align uBAM reads to indexed reference and emit aligned bam output
        java -Xmx"$_MEM_HALF" -jar /usr/gitc/picard.jar SamToFastq \
            INPUT=${reads_bam} \
            FASTQ=/dev/stdout \
            INTERLEAVE=true \
            VALIDATION_STRINGENCY=LENIENT \
        | bwa mem -t `nproc` -R ${read_group} -p ${ref_fasta} - \
        | samtools view -bS - \
        > tempfile1.bam

        # sort, markdup, reorder, and index aligned bam
        java -Xmx"$_MEM" -jar /usr/gitc/picard.jar SortSam \
            I=tempfile1.bam \
            O=tempfile2.bam \
            SO=coordinate
        rm tempfile1.bam
        java -Xmx"$_MEM" -jar /usr/gitc/picard.jar MarkDuplicates \
            I=tempfile2.bam \
            O=tempfile3.bam \
            M=marked_duplicates.metrics
        rm tempfile2.bam
        java -Xmx"$_MEM" -jar /usr/gitc/picard.jar ReorderSam \
            I=tempfile3.bam \
            O=${sample_name}.aligned.bam \
            R=${ref}
        rm tempfile3.bam
        samtools index ${sample_name}.aligned.bam
    }
    output {
        File aligned_bam = "${sample_name}.aligned.bam"
        File aligned_bam_idx = "${sample_name}.aligned.bai"
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "14 GB"
        cpu: 8
        disks: "local-disk 100 LOCAL"
    }
}


# base quality score recalibration
task BQSR {
    File ref_fasta
    File ref_dict
    File ref_index
    File bam 
    File bam_index
    String sample_name
    String output_table_name
    String output_bam_name
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] intervals_to_exclude
    
    command {  
        # build BQSR table 
        java -Xmx7G -jar /usr/gitc/GATK36.jar -T BaseRecalibrator -nt 1 \
            -R ${ref_fasta} -I ${bam} \
            -knownSites ${sep=" -knownSites " known_sites} \
            --excludeIntervals ${sep=" --excludeIntervals " intervals_to_exclude} \
            -o ${output_table_name}

        # install GATK AnalyzeCovariates R dependencies
        R --vanilla << CODE
        install.packages("gplots", repos="http://cran.us.r-project.org")
        install.packages("gsalib", repos="http://cran.us.r-project.org")
        install.packages("reshape", repos="http://cran.us.r-project.org")
        CODE

        # AnalyzeCovariates
        java -jar /usr/gitc/GATK36.jar \
            -T AnalyzeCovariates \
            -R ${ref_fasta} \
            --BQSR ${output_table_name} \
            -plots ${sample_name}.bqsr.pdf

        # clean reads, using bqsr if applicable
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T PrintReads \
            -nt 1 \
            -R ${ref_fasta} \
            -I ${bam} \
            --BQSR ${output_table_name} \
            -o ${output_bam_name}

        # build index
        java -Xmx7G -jar /usr/gitc/picard.jar \
            BuildBamIndex \
            I=${output_bam_name} \
            O=${output_bam_name}.bai
    }
    
    output {
        File out = "${output_bam_name}"
        File out_index = "${output_bam_name}.bai"
        File table = "${output_table_name}"
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB" 
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}

# call snp and indel variants 
task HaplotypeCaller {
    File bam
    File ref
    File ref_dict
    File ref_index
    File bam_index
    File bqsr_table
    File intervals
    String sample_name
    String output_name = "${sample_name}.g.vcf"

    command {
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T HaplotypeCaller \
            -R ${ref} \
            --input_file ${bam} \
            --intervals ${intervals} \
            --BQSR ${bqsr_table} \
            -ERC GVCF \
            --interval_padding 100 \
            -o ${output_name} \
            -variant_index_type LINEAR \
            -variant_index_parameter 128000
    }

    output {
		File vcf = "${output_name}"
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
  }

# merge and genotype vcfs 
task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
	File intervals
    Boolean all_sites
	Array[File] vcf_files
    String gcvf_out = "genotypedGVCFs.vcf"

    command {
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T GenotypeGVCFs \
            -R ${ref} \
            --intervals ${sep=" --intervals " intervals} \
            -o ${gcvf_out} \
            -V ${sep=" -V " vcf_files} \
            ${true="-allSites" false="" all_sites}
    }
    output {
	   File out = gcvf_out
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "13 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}

# variant quality score recalibration
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
task VQSR {
	File ref
    File ref_dict
    File ref_index
    File gvcf 
    File intervals
    String output_filename

    String mode
    Float ts_filter
    Array[String] resources
    Array[String] annotations
	Array[File] resource_files
    Array[File] resource_file_indices
    
    Int max_gaussians
    Int mapping_qual_cap

    String vqsr_file = "${mode}.recal"
    String rscript_file = "${mode}.plots.R"
    String tranches_file = "${mode}.tranches"

	command {
        # build vqsr file
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T VariantRecalibrator \
            -R ${ref} \
            --input ${gvcf} \
            --mode ${mode} \
            --recal_file ${vqsr_file} \
            --tranches_file ${tranches_file} \
            --rscript_file ${rscript_file} \
            --intervals ${sep="--intervals " intervals} \
            --resource:${sep=" --resource:" resources} \
            --use_annotation ${sep=" --use_annotation " annotations} \
            --maxGaussians ${max_gaussians} \
            --MQCapForLogitJitterTransform ${mapping_qual_cap}

        # apply vqsr
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T ApplyRecalibration \
            -R ${ref} \
            --input ${gvcf} \
            --ts_filter_level ${ts_filter} \
            --tranches_file ${tranches_file} \
            --recal_file ${vqsr_file} \
            --mode ${mode} \
            -o ${mode}_vqsr.filtered.vcf
    }

    output {
        File vqsr = vqsr_file
        File rscript = rscript_file
		File tranches = tranches_file
        File out = output_filename
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}

# hard-filter a vcf, if vqsr not available 
# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
task HardFiltration {
    File vcf
    File ref
    File ref_dict
    File ref_index
    String output_filename

    String snp_filter_expr
    String indel_filter_expr 

    command {
        # select snps
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels 
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.g.vcf

        # filter indels
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T CombineVariants \
            -R ${ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ${output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "${output_filename}"
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "13 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}

# annotate variants
# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
    File vcf    
    File ref 
    File Pf3D7_gff
    String output_filename

    String Pf3D7_db_dir = "/opt/snpEff/data/Pf3D7/"

    command {
        # init database
        echo "Pf3D7.genome : Plasmodium_falciparum_3D7" >> /opt/snpEff/snpEff.config
        mkdir -p ${Pf3D7_db_dir}
        mv ${ref} ${Pf3D7_db_dir}/sequences.fa
        mv ${Pf3D7_gff} ${Pf3D7_db_dir}/genes.gff

        # build db 
        java -jar /opt/snpEff/snpEff.jar build -gff3 -v Pf3D7

        # run snpeff 
        java -Xmx7G -jar /opt/snpEff/snpEff.jar \
        	-config /opt/snpEff/snpEff.config \
            -formatEff -no-downstream -no-intergenic \
            -no-upstream -no-utr -noStats \
            -treatAllAsProteinCoding false \
            Pf3D7 ${vcf} > ${output_filename}
    }

    output {
        File out = "${output_filename}"
    }

    runtime { 
        docker: "maxulysse/snpeff:1.3"
        memory: "7 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}
