# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
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

    command {
        java -Xms4G -jar ${gatk_jar} \
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
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    String gatk_jar
    String input_bam
    String output_bam_path
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int? compression_level

    command {
        java ${"-Dsamjdk.compression_level=" + compression_level} -Xms4G -jar ${gatk_jar} \
          ApplyBQSR \
          --create-output-bam-md5 \
          --add-output-sam-program-record \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${output_bam_path}.bam \
          -bqsr ${recalibration_report} \
          --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
          -L ${sep=" -L " sequence_group_interval}
    }
    output {
        File recalibrated_bam = "${output_bam_path}"
        File recalibrated_bam_checksum = "${output_bam_path}.md5"
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    String gatk_jar
    Array[File] input_bqsr_reports
    String output_report_filepath

    command {
        java -Xms3G -jar ${gatk_jar} \
        GatherBQSRReports \
        -I ${sep=' -I ' input_bqsr_reports} \
        -O ${output_report_filepath}
    }
    output {
        File output_bqsr_report = "${output_report_filepath}"
    }
}
