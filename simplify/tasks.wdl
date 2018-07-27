# 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#     - simplified to subset of functionality by dpark/Viral, IDMP, Broad Institute.
#  

task IndexFasta {
    # This task exists, but is probably only useful for large (GB) genomes
    # It's certainly not useful for a 24MB genome
    File   ref_fasta

    String file_basename = basename(ref_fasta, '.fasta')
    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        ln -s ${ref_fasta} ${file_basename}.fasta
        samtools faidx ${file_basename}.fasta
        java -Xmx3G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=${file_basename}.dict
        ls -alF ${file_basename}*
    }
    output {
        File  ref_idx_dict     = "${file_basename}.dict"
        File  ref_idx_fai      = "${file_basename}.fasta.fai"
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB"
        cpu:    2
        disks:  "local-disk 100 SSD"
    }
}

task AlignSortDedupReads {
    File reads_bam     # unaligned reads in BAM, SAM, or CRAM format
    File ref_fasta     # reference fasta

    String file_basename = basename(reads_bam, ".bam")

    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # index fasta
        _REF_FASTA_LOCAL=`basename ${ref_fasta}`
        ln -s ${ref_fasta} $_REF_FASTA_LOCAL
        bwa index $_REF_FASTA_LOCAL
        samtools faidx $_REF_FASTA_LOCAL
        java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict
        cat `basename $_REF_FASTA_LOCAL .fasta`.dict

        # get list of read groups
        samtools view -H ${reads_bam} | grep ^@RG | sed -E -e 's/^@RG.*ID:([^[:space:]]+).*$/\1/' | tee read_group_ids.txt

        # align one RG at a time (because BWA can only write
        # metadata for one RG at a time)
        touch bam_filenames.txt
        for _RG_ID in `cat read_group_ids.txt`; do
            echo "aligning read group $_RG_ID"

            # strip input bam to only one read group
            samtools view -b1 -r $_RG_ID ${reads_bam} > tempfile0.bam
            _RG_DATA=`samtools view -H ${reads_bam} | grep ^@RG | grep ID:$_RG_ID | sed -E -e 's/[[:space:]]+/\\\\t/g'`

            # TO DO: defend against empty tempfile0.bam (zero reads for this RG)

            # align uBAM reads to indexed reference and emit aligned bam output
            java -Xmx2G -jar /usr/gitc/picard.jar SamToFastq \
                INPUT=tempfile0.bam \
                FASTQ=/dev/stdout \
                INTERLEAVE=true \
                VALIDATION_STRINGENCY=LENIENT | \
            bwa mem -t `nproc` -R "$_RG_DATA" -p $_REF_FASTA_LOCAL - | \
            samtools view -b1S - > \
            tempfile1.bam
            rm tempfile0.bam

            # sort bam and record filename in text file
            samtools sort -@ `nproc` -l 1 -T tmpsort -o aligned_$_RG_ID.bam tempfile1.bam
            rm tempfile1.bam
            echo "I=aligned_$_RG_ID.bam" >> bam_filenames.txt
        done

        # TO DO: harden generally against empty input bam (zero reads)

        # merge bwa results for each read group
        java -Xmx12G -jar /usr/gitc/picard.jar MergeSamFiles \
            SORT_ORDER=coordinate USE_THREADING=true CREATE_INDEX=true \
            `cat bam_filenames.txt` \
            O=temp_merged_aligned.bam

        # MarkDuplicates and index
        java -Xmx12G -jar /usr/gitc/picard.jar MarkDuplicates \
            I=temp_merged_aligned.bam \
            O=${file_basename}.aligned.bam \
            M=${file_basename}.aligned.markdup.metrics
        rm temp_merged_aligned.bam
        samtools index ${file_basename}.aligned.bam ${file_basename}.aligned.bai

        # collect figures of merit
        cat `basename $_REF_FASTA_LOCAL .fasta`.dict | grep "^@SQ" | sed -E -e 's/^@SQ.*LN:([^[:space:]]+).*$/\1/' | python -c 'import sys; print(sum(int(l) for l in sys.stdin))' | tee ref_length
        #grep -v '^>' ${ref_fasta}.fasta | tr -d '\n' | wc -c | tee ref_length
        samtools view -c ${file_basename}.aligned.bam | tee reads_aligned
        samtools flagstat ${file_basename}.aligned.bam | tee ${file_basename}.aligned.flagstat.txt
        grep properly ${file_basename}.aligned.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
        samtools view -q 1 -F 1028 ${file_basename}.aligned.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
        python -c "print (float("`cat bases_aligned`")/"`cat ref_length`") if "`cat ref_length`">0 else 0" > mean_coverage
    }
    output {
        File  aligned_bam          = "${file_basename}.aligned.bam"
        File  aligned_bam_idx      = "${file_basename}.aligned.bai"
        File  aligned_bam_flagstat = "${file_basename}.aligned.flagstat.txt"
        File  markdup_metrics      = "${file_basename}.aligned.markdup.metrics"
        Int   reads_aligned        = read_int("reads_aligned")
        Int   read_pairs_aligned   = read_int("read_pairs_aligned")
        Int   bases_aligned        = read_int("bases_aligned")
        Int   ref_length           = read_int("ref_length")
        Float mean_coverage        = read_float("mean_coverage")
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "14 GB"
        cpu:    8
        disks:  "local-disk 100 SSD"
    }
}


# base quality score recalibration
task BaseRecalibrator {
    File           ref_fasta
    File           aligned_bam
    File           aligned_bam_idx

    Array[File]?   known_sites
    Array[String]? intervals_to_exclude

    String         file_basename = basename(aligned_bam, ".bam")
    
    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # prep known_sites/exclude params for BQSR
        BQSR_KNOWN_SITES="${sep=' -knownSites ' known_sites}"
        BQSR_EXCLUDE_INTERVALS="${sep=' --excludeIntervals ' intervals_to_exclude}"
        if [ -n "$BQSR_KNOWN_SITES" ]; then BQSR_KNOWN_SITES="-knownSites $BQSR_KNOWN_SITES"; fi
        if [ -n "$BQSR_EXCLUDE_INTERVALS" ]; then BQSR_EXCLUDE_INTERVALS="--excludeIntervals $BQSR_EXCLUDE_INTERVALS"; fi

        if [ -n "$BQSR_KNOWN_SITES" ]; then
            # known_sites is supplied, so run BQSR

            # index fasta
            _REF_FASTA_LOCAL=`basename ${ref_fasta}`
            ln -s ${ref_fasta} $_REF_FASTA_LOCAL
            samtools faidx $_REF_FASTA_LOCAL
            java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
                R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict

            # build BQSR table 
            java -Xmx7G -jar /usr/gitc/GATK36.jar -T BaseRecalibrator \
                -nct `nproc` \
                -R $_REF_FASTA_LOCAL -I ${aligned_bam} \
                $BQSR_KNOWN_SITES $BQSR_EXCLUDE_INTERVALS \
                -o ${file_basename}.bqsr.txt

            # rewrite qual scores to new BAM
            java -Xmx3G -jar /usr/gitc/GATK36.jar \
                -T PrintReads \
                -nct `nproc` \
                -R $_REF_FASTA_LOCAL \
                -I ${aligned_bam} \
                --BQSR ${file_basename}.bqsr.txt \
                -o ${file_basename}.bqsr.bam &

            # AnalyzeCovariates (depends on R libraries)
            echo 'install.packages("gplots", repos="http://cran.us.r-project.org"); install.packages("gsalib", repos="http://cran.us.r-project.org"); install.packages("reshape", repos="http://cran.us.r-project.org")' | R --vanilla
            java -jar -Xmx3G /usr/gitc/GATK36.jar \
                -T AnalyzeCovariates \
                -R $_REF_FASTA_LOCAL \
                --BQSR ${file_basename}.bqsr.txt \
                -plots ${file_basename}.bqsr.pdf

            wait # for PrintReads to finish
            # index recalibrated aligned bam
            samtools index ${file_basename}.bqsr.bam ${file_basename}.bqsr.bai

        else
            # known_sites is empty, therefore skip BQSR (which requires known_sites)
            # cp instead of ln -s because of how cromwell handles chroot/docker stuff
            cp ${aligned_bam} ${file_basename}.bqsr.bam
            cp ${aligned_bam_idx} ${file_basename}.bqsr.bai
            touch ${file_basename}.bqsr.txt ${file_basename}.bqsr.pdf

        fi
    }
    
    output {
        File table       = "${file_basename}.bqsr.txt"
        File out_bam     = "${file_basename}.bqsr.bam"
        File out_bam_idx = "${file_basename}.bqsr.bai"
        File bqsr_plot   = "${file_basename}.bqsr.pdf"
    } 

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        memory: "7 GB" 
        cpu:    2
        disks:  "local-disk 100 SSD"
    }
}


# call snp and indel variants 
task HaplotypeCaller {
    File   aligned_bam
    File   aligned_bam_idx
    File   ref_fasta

    File?  bqsr_table
    File?  intervals   # turn this into Array[String] perhaps?
    Int?   ploidy
    Int?   min_base_quality_score
    Int?   downsample_to_coverage
    Int?   interval_padding

    String file_basename = basename(aligned_bam, ".bam")

    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # index fasta
        _REF_FASTA_LOCAL=`basename ${ref_fasta}`
        ln -s ${ref_fasta} $_REF_FASTA_LOCAL
        samtools faidx $_REF_FASTA_LOCAL
        java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict

        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T HaplotypeCaller \
            -R $_REF_FASTA_LOCAL \
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
        cpu:    2
        disks:  "local-disk 100 SSD"
    }
  }

# merge and genotype vcfs 
task GenotypeGVCFs {
    Array[File]+ vcf_files
    Array[File]+ vcf_indexes_tbi
    File         ref_fasta

    File?        intervals  # convert to Array[String] perhaps?
    Boolean      all_sites=true
    String       out_basename = "joint_genotypes"

    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # index fasta
        _REF_FASTA_LOCAL=`basename ${ref_fasta}`
        ln -s ${ref_fasta} $_REF_FASTA_LOCAL
        samtools faidx $_REF_FASTA_LOCAL
        java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict

        java -Xmx12G -jar /usr/gitc/GATK36.jar \
            -T GenotypeGVCFs \
            -R $_REF_FASTA_LOCAL \
            -V ${sep=" -V " vcf_files} \
            ${'--intervals ' + intervals} \
            ${true="-allSites" false="" all_sites} \
            -nt `nproc` \
            --out ${out_basename}.vcf
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
        disks: "local-disk 100 SSD"
    }
}

# variant quality score recalibration
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
task VQSR {
    File          ref_fasta
    File          gvcf 
    File          intervals
    String        output_filename

    String        mode
    Float         ts_filter
    Array[String] resources
    Array[String] annotations
    Array[File]   resource_files
    Array[File]   resource_file_indices
    
    Int max_gaussians
    Int mapping_qual_cap

    # vqsr default values for SNP vs INDEL
    #Float ts_filter_snp=99.0
    #Float ts_filter_indel=99.5
    #Int snp_max_gaussians=8
    #Int indel_max_gaussians=4
    #Int vqsr_mapping_qual_cap=70
    #Array[String] snp_resources=["7g8_gb4,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=15.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/3d7_hb3.combined.final.vcf.gz"]
    #Array[String] indel_resources=["7g8_gb4,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=12.0 /cromwell_root/fc-6d3db020-0f34-4434-a132-b1a4719a6848/refs/3d7_hb3.combined.final.vcf.gz"]
    #Array[String] snp_annotations=["QD", "FS", "SOR", "DP"]
    #Array[String] indel_annotations=["QD", "FS"]

    String vqsr_file = "${mode}.recal"
    String rscript_file = "${mode}.plots.R"
    String tranches_file = "${mode}.tranches"

    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # index fasta
        _REF_FASTA_LOCAL=`basename ${ref_fasta}`
        ln -s ${ref_fasta} $_REF_FASTA_LOCAL
        samtools faidx $_REF_FASTA_LOCAL
        java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict

        # build vqsr file
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T VariantRecalibrator \
            -R $_REF_FASTA_LOCAL \
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
            -R $_REF_FASTA_LOCAL \
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
        disks: "local-disk 100 SSD"
    }
}

# hard-filter a vcf, if vqsr not available 
# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
task HardFiltration {
    File   vcf
    File?  vcf_tbi
    File   ref_fasta

    String snp_filter_expr = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    String indel_filter_expr = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

    String base_filename = basename(basename(vcf, ".gz"), ".vcf")

    command {
        set -ex -o pipefail
        PATH=$PATH:/usr/gitc

        # index fasta
        _REF_FASTA_LOCAL=`basename ${ref_fasta}`
        ln -s ${ref_fasta} $_REF_FASTA_LOCAL
        samtools faidx $_REF_FASTA_LOCAL
        java -Xmx2G -jar /usr/gitc/picard.jar CreateSequenceDictionary \
            R=$_REF_FASTA_LOCAL O=`basename $_REF_FASTA_LOCAL .fasta`.dict

        # select & filter snps -> filtered_snps.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R $_REF_FASTA_LOCAL \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R $_REF_FASTA_LOCAL \
            -V raw_snps.vcf \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.vcf
        rm raw_snps.vcf

        # select & filter indels -> filtered_indels.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T SelectVariants \
            -R $_REF_FASTA_LOCAL \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.vcf
        java -Xmx3G -jar /usr/gitc/GATK36.jar \
            -T VariantFiltration \
            -R $_REF_FASTA_LOCAL \
            -V raw_indels.vcf \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.vcf
        rm raw_indels.vcf

        # combine variants
        java -Xmx7G -jar /usr/gitc/GATK36.jar \
            -T CombineVariants \
            -R $_REF_FASTA_LOCAL \
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
        disks: "local-disk 100 SSD"
    }
}

# annotate variants
# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
    File   vcf
    File?  vcf_tbi
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
        disks: "local-disk 100 SSD"
    }
}
