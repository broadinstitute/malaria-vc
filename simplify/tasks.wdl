# 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

## TASK DEFINITIONS 

task AlignSortDedupReads {
    File reads_bam     # unaligned reads in BAM, SAM, or CRAM format
    File ref_fasta     # not 100% sure this input file is necessary
    File ref_idx_dict  # Picard/GATK index
    File ref_idx_fai   # samtools index
    File ref_idx_amb   # bwa index
    File ref_idx_ann   # bwa index
    File ref_idx_bwt   # bwa index
    File ref_idx_pac   # bwa index
    File ref_idx_sa    # bwa index

    String file_basename = basename(reads_bam, ".bam")

    command {
        set -ex -o pipefail
        _MEM="12G"
        _MEM_HALF="6G"

        # TODO: pull {read_group} from in_bam header

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
            O=${file_basename}.aligned.bam \
            R=${ref}
        rm tempfile3.bam
        samtools index ${file_basename}.aligned.bam
    }
    output {
        File aligned_bam = "${file_basename}.aligned.bam"
        File aligned_bam_idx = "${file_basename}.aligned.bai"
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
    File ref_idx_dict
    File ref_idx_fai
    File aligned_bam 
    File aligned_bam_idx
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] intervals_to_exclude

    String file_basename = basename(aligned_bam, ".bam")
    
    command {
        set -ex -o pipefail

        # build BQSR table 
        java -Xmx7G -jar /usr/gitc/GATK36.jar -T BaseRecalibrator -nt 1 \
            -R ${ref_fasta} -I ${aligned_bam} \
            -knownSites ${sep=" -knownSites " known_sites} \
            --excludeIntervals ${sep=" --excludeIntervals " intervals_to_exclude} \
            -o ${file_basename}.bqsr.txt

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
            --BQSR ${file_basename}.bqsr.txt \
            -plots ${file_basename}.bqsr.pdf

        # clean reads, using bqsr if applicable
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T PrintReads \
            -nt 1 \
            -R ${ref_fasta} \
            -I ${aligned_bam} \
            --BQSR ${file_basename}.bqsr.txt \
            -o ${file_basename}.bqsr.bam

        # build index
        java -Xmx7G -jar /usr/gitc/picard.jar \
            BuildBamIndex \
            I=${file_basename}.bqsr.bam \
            O=${file_basename}.bqsr.bai
    }
    
    output {
        File out_bam = "${file_basename}.bqsr.bam"
        File out_bam_idx = "${file_basename}.bqsr.bai"
        File table = "${file_basename}.bqsr.txt"
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
    File   aligned_bam
    File   aligned_bam_idx
    File   ref_fasta
    File   ref_idx_dict
    File   ref_idx_fai

    File?  bqsr_table
    File?  intervals   # turn this into Array[String] perhaps?
    Int?   ploidy
    Int?   min_base_quality_score
    Int?   downsample_to_coverage
    Int?   interval_padding

    String file_basename = basename(aligned_bam, ".bam")

    command {
        set -ex -o pipefail
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T HaplotypeCaller \
            -R ${ref_fasta} \
            --input_file ${aligned_bam} \
            --emitRefConfidence GVCF \
            ${'--intervals ' + intervals} \
            ${'--BQSR ' + bqsr_table} \
            ${'--interval_padding ' + interval_padding} \
            ${'--sample_ploidy ' + ploidy} \
            ${'--downsample_to_coverage ' + downsample_to_coverage} \
            ${'--min_base_quality_score ' + min_base_quality_score} \
            --out ${file_basename}.g.vcf
        /usr/gitc/bgzip -c --threads `nproc` ${file_basename}.g.vcf > ${file_basename}.g.vcf.gz
        /usr/gitc/tabix -p vcf ${file_basename}.g.vcf.gz
    }

    output {
        File vcf     = "${file_basename}.g.vcf.gz"
        File vcf_tbi = "${file_basename}.g.vcf.gz.tbi"
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
    Array[File]+ vcf_files
    File         ref_fasta
    File         ref_idx_dict
    File         ref_idx_fai

    File?        intervals  # convert to Array[String] perhaps?
    Boolean      all_sites=true
    String       out_basename = "joint_genotypes"

    command {
        set -ex -o pipefail
        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T GenotypeGVCFs \
            -R ${ref_fasta} \
            -V ${sep=" -V " vcf_files} \
            ${'--intervals ' + intervals} \
            ${true="-allSites" false="" all_sites} \
            -nt `nproc` \
            --out ${gcvf_out}
        /usr/gitc/bgzip -c --threads `nproc` ${out_basename}.vcf > ${out_basename}.vcf.gz
        /usr/gitc/tabix -p vcf ${out_basename}.vcf.gz
    }
    output {
       File out_vcf     = "${out_basename}.vcf.gz"
       File out_vcf_tbi = "${out_basename}.vcf.gz.tbi"
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
    File         ref_fasta
    File         ref_idx_dict
    File         ref_idx_fai
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
            -R ${ref_fasta} \
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
            -R ${ref_fasta} \
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
    File   vcf
    File   ref_fasta
    File   ref_idx_dict
    File   ref_idx_fai

    String snp_filter_expr = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    String indel_filter_expr = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

    String base_filename = basename(basename(vcf, ".gz"), ".vcf")

    command {
        set -ex -o pipefail

        # select & filter snps -> filtered_snps.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R ${ref_fasta} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R ${ref_fasta} \
            -V raw_snps.vcf \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.vcf
        rm raw_snps.vcf

        # select & filter indels -> filtered_indels.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R ${ref_fasta} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R ${ref_fasta} \
            -V raw_indels.vcf \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.vcf
        rm raw_indels.vcf

        # combine variants
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T CombineVariants \
            -R ${ref_fasta} \
            --variant filtered_snps.vcf \
            --variant filtered_indels.vcf \
            -o ${base_filename}.filtered.vcf \
            --genotypemergeoption REQUIRE_UNIQUE
        /usr/gitc/bgzip -c --threads `nproc` ${base_filename}.filtered.vcf > ${base_filename}.filtered.vcf.gz
        /usr/gitc/tabix -p vcf ${base_filename}.filtered.vcf.gz
    }

    output {
        File out_vcf     = "${base_filename}.filtered.vcf.gz"
        File out_vcf_tbi = "${base_filename}.filtered.vcf.gz.tbi"
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}

# annotate variants
# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
    File   vcf    
    File   ref_fasta 
    File   ref_gff

    String base_filename = basename(basename(vcf, ".gz"), ".vcf")

    command {
        set -ex -o pipefail

        # init database
        echo "custom_genome.genome : Plasmodium_falciparum_3D7" >> /opt/snpEff/snpEff.config
        mkdir -p /opt/snpEff/data/custom_genome
        mv ${ref_fasta} /opt/snpEff/data/custom_genome/sequences.fa
        mv ${ref_gff} /opt/snpEff/data/custom_genome/genes.gff

        # build db 
        java -jar /opt/snpEff/snpEff.jar build -gff3 -v custom_genome

        # run snpeff 
        java -Xmx7G -jar /opt/snpEff/snpEff.jar \
            -config /opt/snpEff/snpEff.config \
            -formatEff -no-downstream -no-intergenic \
            -no-upstream -no-utr -noStats \
            -treatAllAsProteinCoding false \
            custom_genome \
            ${vcf} \
            > ${base_filename}.snpeff.vcf
    }

    output {
        File annotated_vcf = "${base_filename}.snpeff.vcf"
    }

    runtime { 
        docker: "maxulysse/snpeff:1.3"
        memory: "7 GB"
        cpu: 2
        disks: "local-disk 100 LOCAL"
    }
}
