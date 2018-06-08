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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    String? preCommand
    Array[File]+ inputBams
    Array[File]+ inputBamsIndex
    Array[File]+ intervalList
    String gvcfPath
    File refDict
    File refFasta
    File refFastaIndex
    Float? contamination
    Int? compressionLevel
    String gatkJar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compressionLevel} \
        -Xmx${mem}G -jar ${gatkJar} \
          HaplotypeCaller \
          -R ${refFasta} \
          -O ${gvcfPath} \
          -I ${sep=" -I " inputBams} \
          -L ${sep=' -L ' intervalList} \
          -contamination ${default=0 contamination} \
          -ERC GVCF
    }

    output {
        File outputGVCF = gvcfPath
        File outputGVCFindex = gvcfPath + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task GenotypeGVCFs {
    String? preCommand
    File gvcfFiles
    File gvcfFileIndexes
    Array[File]+ intervals

    String outputPath

    String gatkJar

    File refFasta
    File refFastaIndex
    File refDict

    File dbsnpVCF
    File dbsnpVCFindex

    Int? compressionLevel
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        java ${"-Dsamjdk.compression_level=" + compressionLevel} \
        -Xmx${mem}G -jar ${gatkJar} \
         GenotypeGVCFs \
         -R ${refFasta} \
         -O ${outputPath} \
         -D ${dbsnpVCF} \
         -G StandardAnnotation \
         --only-output-calls-starting-in-intervals \
         -new-qual \
         -V ${gvcfFiles} \
         -L ${sep=' -L ' intervals}
    }

    output {
        File outputVCF = outputPath
        File outputVCFindex = outputPath + ".tbi"
    }

    runtime{
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task CombineGVCFs {
    String? preCommand
    Array[File]+ gvcfFiles
    Array[File]+ gvcfFileIndexes
    Array[File]+ intervals

    String outputPath

    String gatkJar

    File refFasta
    File refFastaIndex
    File refDict

    Int? compressionLevel
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        if [ ${length(gvcfFiles)} -gt 1 ]; then
            java ${"-Dsamjdk.compression_level=" + compressionLevel} \
            -Xmx${mem}G -jar ${gatkJar} \
             CombineGVCFs \
             -R ${refFasta} \
             -O ${outputPath} \
             -V ${sep=' -V ' gvcfFiles} \
             -L ${sep=' -L ' intervals}
        else
            ln -sf ${select_first(gvcfFiles)} ${outputPath}
            ln -sf ${select_first(gvcfFileIndexes)} ${outputPath}.tbi
        fi
    }

    output {
        File outputGVCF = outputPath
        File outputGVCFindex = outputPath + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task SplitNCigarReads {
    String? preCommand

    File inputBam
    File inputBamIndex
    File refFasta
    File refFastaIndex
    File refDict
    String outputBam
    String gatkJar
    Array[File]+ intervals

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatkJar} \
        SplitNCigarReads \
        -I ${inputBam} \
        -R ${refFasta} \
        -O ${outputBam} \
        -L ${sep=' -L ' intervals}
    }

    output {
        File bam = outputBam
        File bam_index = sub(outputBam, "\\.bam$", ".bai")
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}
