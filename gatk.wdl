# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    String? preCommand
    String gatk_jar
    String input_bam
    String input_bam_index
    String recalibration_report_filename
    Array[File]+ sequence_group_interval
    Array[File]+ known_indels_sites_VCFs
    Array[File]+ known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatk_jar} \
          BaseRecalibrator \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${recalibration_report_filename} \
          --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
          -L ${sep=" -L " sequence_group_interval}
    }

    output {
        File recalibration_report = "${recalibration_report_filename}"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    String? preCommand
    String gatk_jar
    String input_bam
    String output_bam_path
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int? compression_level

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xms${mem}G -jar ${gatk_jar} \
          ApplyBQSR \
          --create-output-bam-md5 \
          --add-output-sam-program-record \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${output_bam_path} \
          -bqsr ${recalibration_report} \
          --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
          -L ${sep=" -L " sequence_group_interval}
    }

    output {
        File recalibrated_bam = "${output_bam_path}"
        File recalibrated_bam_checksum = "${output_bam_path}.md5"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    String? preCommand
    String gatk_jar
    Array[File] input_bqsr_reports
    String output_report_filepath

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatk_jar} \
        GatherBQSRReports \
        -I ${sep=' -I ' input_bqsr_reports} \
        -O ${output_report_filepath}
    }

    output {
        File output_bqsr_report = "${output_report_filepath}"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    String? preCommand
    Array[File]+ input_bams
    Array[File]+ input_bams_index
    Array[File]+ interval_list
    String gvcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Int? compression_level
    String gatk_jar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xmx${mem}G -jar ${gatk_jar} \
          HaplotypeCaller \
          -R ${ref_fasta} \
          -O ${gvcf_basename}.vcf.gz \
          -I ${sep=" -I " input_bams} \
          -L ${sep=' -L ' interval_list} \
          -contamination ${default=0 contamination} \
          -ERC GVCF
    }

    output {
        File output_gvcf = "${gvcf_basename}.vcf.gz"
        File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

task GenotypeGVCFs {
    String? preCommand
    File gvcf_files
    File gvcf_file_indexes
    Array[File]+ intervals

    String output_basename

    String gatk_jar

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int? compression_level
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xmx${mem}G -jar ${gatk_jar} \
         GenotypeGVCFs \
         -R ${ref_fasta} \
         -O ${output_basename + ".vcf.gz"} \
         -D ${dbsnp_vcf} \
         -G StandardAnnotation \
         --only-output-calls-starting-in-intervals \
         -new-qual \
         -V ${gvcf_files} \
         -L ${sep=' -L ' intervals}
    }

    output {
        File output_vcf = output_basename + ".vcf.gz"
        File output_vcf_index = output_basename + ".vcf.gz.tbi"
    }

    runtime{
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

task CombineGVCFs {
    String? preCommand
    Array[File]+ gvcf_files
    Array[File]+ gvcf_file_indexes
    Array[File]+ intervals

    String output_basename

    String gatk_jar

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int? compression_level
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        if [ ${length(gvcf_files)} -gt 1 ]; then
            java ${"-Dsamjdk.compression_level=" + compression_level} \
            -Xmx${mem}G -jar ${gatk_jar} \
             CombineGVCFs \
             -R ${ref_fasta} \
             -O ${output_basename + ".vcf.gz"} \
             -V ${sep=' -V ' gvcf_files} \
             -L ${sep=' -L ' intervals}
        else
            ln -sf ${select_first(gvcf_files)} ${output_basename + ".vcf.gz"}
            ln -sf ${select_first(gvcf_files)}.tbi ${output_basename + ".vcf.gz.tbi"}
        fi
    }

    output {
        File output_gvcf = output_basename + ".vcf.gz"
        File output_gvcf_index = output_basename + ".vcf.gz.tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

task SplitNCigarReads {
    String? preCommand

    File input_bam
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String output_bam
    String gatk_jar
    Array[File]+ intervals

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatk_jar} \
        SplitNCigarReads \
        -I ${input_bam} \
        -R ${ref_fasta} \
        -O ${output_bam} # might have to be -o depending on GATK version \
        -L ${sep=' -L ' intervals}
    }

    output {
        File bam = output_bam
        File bam_index = output_bam + ".bai"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}
